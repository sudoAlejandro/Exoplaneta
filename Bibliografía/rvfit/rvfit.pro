;+
; NAME:
;	rvfit
;
; PURPOSE:
;	Fits the parameters for a keplerian radial velocity function using ASA.
;
; EXPLANATION:
;	rvfit fits a keplerian radial velocity function to get the values for the parameters
;	P, Tp, e, omega, gamma, K1, K2. If the system is double-lined all the parameters are
;	fitted. If the system is single-lined K2 is not fitted. This code uses ASA, Adaptive
;	Simulated Annealing to fit simultaneously all the parameters.
;
; CALLING SEQUENCE:
;	rvfit,configfile=configfile[,outfile=outfile][,evolfile=evolfile][,/physics][,/autodom]
;	rvfit,rvfile1,rvfile2,fitparam,valparam,L,U[,outfile=outfile][,evolfile=evolfile][,/physics][,/autodom]
;
; INPUTS:
;	configfile = name of the configuration file (string). It must contain the values for
;		the variables rvfile1, rvfile2, fitparam, valparam, LL, and UU, one per line.
;	outputfile = name of the output file (string).
;	evolfile = name of the evolution file (string). It contains the parameter values as
;		the algorithm evolves in the parameter space.
;	rvfile1 = name of the file containing the primary radial velocities (string).
;	rvfile2 = name of the file containing the secondary radial velocities (string).
;	fitparam = seven elements vector with 0 or 1 depending wether the param is fitted.
;		The parameters must be ordered following [P, Tp, e, omega, gamma, K1, K2].
;	valparam = seven elements vector with the values for the non fitted parameters, with
;		a 0 in fitparam.
;	L = seven elements vector with the lower limit of the paramters.
;	U = seven elements vector with the upper limit of the paramters.
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;	/physics = compute the physical parameters of the system and propagate their
;		uncertainties.
;	/autodom = compute the domain limits, i.e. the L and U vectors, based in the data.
;
; OUTPUTS
;	Fitted values for the parameters [P, Tp, e, omega, gamma, K1, K2] and their uncertainties.
;	Plot with the RV curve and the residuals. Values for the physical parameters.
;
; EXAMPLES
;	To fit all the keplerian parameters of LV Her and compute the physical parameters:
;	rvfit,'LVHer_RV1.dat','LVHer_RV2.dat', $
;		[1,1,1,1,1,1,1], $
;		[18.4359535,2453652.19147d,0.61273,352.20,-10.278,67.24d,68.59d], $
;		[10.,2448407.8794d,0.,0.,-100.,0.,0.], $
;		[ 25.,2452425.8862d,0.999,360.,100.,100.,150.],/physics
;
;	It is easier to put all the parameters in a config file and run rvfit with it. In this
;	case we also want to see how the parameters evolves:
;	rvfit,configfile='LVHer.conf',outfile='LVHer_param.out',evolfile='LVHer_ASA.out',/physics
;
;	A typical config file would contain this lines (whithout the ;):
;
; rvfile1  = LVHer_RV1.dat
; rvfile2  = LVHer_RV2.dat
; fitparam = [0,          1,              1,       1,       1,      1,      1]
; valparam = [18.4359535, 2448414.90436d, 0.61273, 352.20, -10.278, 67.24d, 68.59d]
; L        = [10.,        2448407.8794d,  0.,      0.,     -100.,   0.,     0.]
; U        = [25.,        2452425.8862d,  0.999,   360.,    100,    100.,   150.]
;
; REFERENCES:
;	S. Chen, B.L. Luk, "Adaptive simulated annealing for optimization in
;		signal processing applications", Signal Processing, 1999, 79, p117-128
;	L. Ingber, "Very fast simulated re-annealing", Mathl. Comput. Modelling, 1989,
;		12, p967-973.
;	L. Ingber, "Simulated annealing: Practice versus theory", Mathl. Comput. Modelling,
;		1993, 18, Nº11, p29-57.
;	L. Ingber, "Adaptive simulated annealing (ASA): Lessons learned", Control and
;		Cybernetics, 1996, Vol 25, Nº1, p33-54.
;	Dreo, Petrowski, Siarry, Taillard, "Metaheuristics For Hard Optimization", Springer.
;
; MODIFICATION HISTORY:
;       Written by Ramon Iglesias Marzoa, Dep. Astrofisica (ULL), 2013-03-09.
;	Last change: 2014-12-30.
;-

; Append the necesary functions and procedures.
@RVlib
@computeRVparams

; Generate a new vector xnew from x using the generation temperatures.
function genera,x,fitparam,Tgen
common limits,L,U
nx=n_elements(x)
xnew=x

for g=0,nx-1 do begin	; This loop sweep across the dimensions.

  ; If this parameter is fixed then continue with the next one.
  if fitparam[g] eq 1 then begin

    ; Xenera novos valores da función dentro do dominio.
    repeat begin
      rand=randomu(seed1,/double)	; rand between 0 and 1.
      sgn=(rand-0.5) ge 0. ? 1. : -1.	; sign of rand-1/2
      q=sgn*Tgen[g]*(( 1.0d + 1.0d /Tgen[g])^(abs(2*rand-1.))-1.)
      xnew[g]=x[g]+q*(U[g]-L[g])
    endrep until xnew[g] ge L[g] and xnew[g] le U[g]
  ;       print,'g=',g,' rand=',rand,' sgn=',sgn,' q=',q, $
  ;       	' U[g]=',U[g],' L[g]=',L[g],' xnew[g]=',xnew[g]
  
  endif
endfor	; End of the loop of tries in each dimension.

; Normalize Tp in the L[1] - L[1]+P interval by doing:
; ciclo=floor((x-x1)/P), xnorm=x-ciclo*P
xnew[1]=xnew[1]-floor((xnew[1]-L[1])/xnew[0])*xnew[0]

return,xnew
end

pro rvfit,rvfile1,rvfile2,fitparam,valparam,LL,UU, $
  configfile=configfile, $
  outfile=outfile, $
  autodom=autodom, $
  evolfile=evolfile, $
  physics=physics

; To measure the execution time.
tprocini = SYSTIME(1)

; Global variables.
common datos,t1,rv1,srv1,t2,rv2,srv2
common limits,L,U

usageconfigfile='rvfit,configfile=configfile[,outfile=outfile][,evolfile=evolfile][,/physics][,/autodom]'
usageparameters='rvfit,rvfile1,rvfile2,fitparam,valparam,L,U[,outfile=outfile][,evolfile=evolfile][,/physics][,/autodom]'

; Distinguish between the two input modes and format the input.
case n_params() of
  0: begin

    ; Check the configfile.
    if ~keyword_set(configfile) then begin
      print,'You must use a valid configuration file:'
      print,usageconfigfile
      retall
    endif
    
    ; Open the configfile and data formatting.
    data=strarr(1,6)	; Number of lines in the configfile.
    openr,uni,configfile,/GET_LUN
    readf,uni,data
    close,uni
    FREE_LUN,uni

    rvfile1=strtrim(strmid(data[0],strpos(data[0],'=')+1),2)
    rvfile2=strtrim(strmid(data[1],strpos(data[1],'=')+1),2)
    temp=strsplit(data[2],'=',/extract)
    fitparam=fix(strsplit(temp[1],' [,]',/extract))
    temp=strsplit(data[3],'=',/extract)
    valparam=double(strsplit(temp[1],' [,]',/extract))
    temp=strsplit(data[4],'=',/extract)
    L=double(strsplit(temp[1],' [,]',/extract))
    temp=strsplit(data[5],'=',/extract)
    U=double(strsplit(temp[1],' [,]',/extract))
  end
  6: begin
    L=LL
    U=UU
  end
  else: begin
    print,'Usage:'
    print,usageconfigfile
    print,usageparameters
    retall
  end
endcase

; Name of the parameters file.
if ~keyword_set(outfile) then outfile='rvfit.out'

; Flag to distinguish between double and single-line binaries.
flagDL=rvfile2 eq '' ? 0 : 1
npar=7	; total number of parameters.

; Check the number of elements.
if n_elements(fitparam) ne npar $
  or n_elements(valparam) ne npar $
  or n_elements(U) ne npar $
  or n_elements(L) ne npar then begin
    print,'fitparam, valparam, L and U must have '+strtrim(npar,2)+' elements.'
    print,'Nothing is done.'
    retall
endif

; Check the domain.
ind=where(U-L le 0,count)
if count gt 0 then begin
  print,'All the elements of U must be greater than of L.'
  print,'Nothing is done.'
  retall
endif

; Check the fitting flags.
ind=where(fitparam eq 1,count)
if count lt 2 then begin
  print,'The number of fitted parameters must be greater or equal to 2.'
  print,'Nothing is done.'
  retall
endif

; Data reading.
readcol,rvfile1,t1,rv1,srv1,f='d,f,f',comment='#',/silent
nobs1=n_elements(rv1)

; Check the uncertainties.
kk=where(srv1 le 0.,count)
if count gt 0 then begin
  print,'All the elements of srv1 must be greater than 0.'
  print,'Nothing is done.'
  retall
endif

; For double-line binaries.
if flagDL eq 1 then begin

  ; Read the secondary RVs.
  readcol,rvfile2,t2,rv2,srv2,f='d,f,f',comment='#',/silent
  nobs2=n_elements(rv2)

  ; Check the uncertainties.
  kk=where(srv2 le 0.,count)
  if count gt 0 then begin
    print,'All the elements of srv2 must be greater than 0.'
    print,'Nothing is done.'
    retall
  endif
  
  tt=[t1,t2]
  rrvv=[rv1,rv2]
  ssrrvv=[srv1,srv2]

  ; Name of the function for single-line binaries.
  funcname='chisqDL'
  
;   chislimit=total(srv1)+total(srv2)
  
  endif else begin
  
  ; For single-line binaries and exoplanets.
  tt=t1
  rrvv=rv1
  ssrrvv=srv1
  
  ; This ensures that K2 is fixed to 0.
  fitparam[6]=0
  valparam[6]=0
  L[6]=0.	; This is to maintain coherence.
  U[6]=1.
  nobs2=0
  
  ; Name of the function for single-line binaries.
  funcname='chisqSL'
  
;  chislimit=total(srv1)
  
endelse

; This limit is to stop the annealing loop.
meansrv=mean(ssrrvv)
chislimit=total((meansrv/ssrrvv)^2)
; print,chislimit
; print,nobs1+nobs2

; If e=0.0 and fixed then set omega=90. and fixed. 
if fitparam[2] eq 0 and valparam[2] eq 0. then begin
  fitparam[3]=0
  valparam[3]=90.
endif

; ------------------Function domain------------------

if keyword_set(autodom) then begin

  ; Compute maximum and minimum period.
  tprov=tt(sort(tt))
  dt=tprov-shift(tprov,1)	; differences among neighbour times.
  Pmin=2.*min(abs(dt))		; Pmin is the inverse of the Nyquist frecuency.
  Pmax=(max(tt)-min(tt))/2.
  
  ;  P,		Tp,		e,	omega,	gamma,		K1,			K2				 
  L=[Pmin,	min(tt),	0.,	0.,	min(rrvv),	0.,			0.]	; lower limit.
  U=[Pmax,	min(tt)+Pmax,	0.999,	360.,	max(rrvv),	max(rv1)-min(rv1),	1.]	; upper limit.
  
  ; Update the K2 limits in case of double-line binary.
  if flagDL eq 1 then begin
    L[6]=0.
    U[6]=max(rv2)-min(rv2)
  endif
  
endif

; ------------------Initial parameter values-----------------

; Starting parameters.
x=(U+L)/2.
;print,transpose([[L],[U],[x]])

; ------------------Initialization---------------------------

; Freeze the parameters set as fixed.
indfixed=where(fitparam eq 0,nfixed,complement=indfitted,ncomplement=nfitted)
if nfixed gt 0 then x[indfixed]=valparam[indfixed]

; First values for f, xbest and fbest.
f=call_function(funcname,x)
xbest=x
fbest=f

; Stopping parameters.
eps=1e-5	; Allowed tolerance in the function minimum value.
Neps=5		; Number of times that tolerance eps is achieved before termination.
Nterm=20  ;npar	; Number of consecutive re-annealings to stop.
fbestlist=replicate(1.,Neps)	; List with the Neps last values of fbest.
nrean=0		; Re-annealing counter.

; Acceptance temperature, it depends on ka.
ka=0UL
Ta0=f
Ta=Ta0
nacep=0UL	; Acceptance counter.

; Generating temperature for each adjusted parameter,
; it depends on kgen (also for each parameter).
kgen=replicate(0UL,npar)
Tgen0=replicate(1.0d,npar)
Tgen=Tgen0

; Adjustable parameters of the algorithm.
; WARNING: changing these parameters can make the algorithm
; doesn't work or becomes too slow.
c=20.
Na=1000L
Ngen=10000L
delta=abs(U-L)*1e-8	; to compute the sensibilities.

; ------------------Initial acceptance temperature-----------------------

; This follows the prescription of Dreo_Petrowski_Siarry_Taillard
; "Metaheuristics For Hard Optimization-Springer", p44.
print,'Setting the initial Ta...'
acep=0.25
ntest=100
ftest=fltarr(ntest)
for j=0L,ntest-1 do begin
  xnew=genera(x,fitparam,Tgen)
  ;print,xnew
  ftest[j]=call_function(funcname,xnew)

  ; Save the best values.
  if ftest[j] lt fbest then begin
    fbest=ftest[j]
    xbest=xnew
;     print,fbest
  endif
endfor
dftest=shift(ftest,-1)-ftest
dftest=dftest[0:ntest-2]
avdftest=mean(abs(dftest))
Ta0=avdftest/alog(1./acep-1.)
print,'Initial Ta = '+strtrim(Ta0,2)

; ------------------Simulated annealing algorithm-----------------------

; Save the Markov Chain in a file.
if keyword_set(evolfile) then begin
  ; Delete previous results.
  ;spawn,'rm -f '+evolfile
  openw,uni,evolfile,/GET_LUN,width=160
  printf,uni,'# Ta F(X) X(7 elements)'
endif

; Annealing loop. It reduces the temperature.
repeat begin

  for j=0L,Ngen-1 do begin	; Loop of generated points for each temperature.

    flag_aceptancia=0

    ; Generate a new value, xnew.
    xnew=genera(x,fitparam,Tgen)
    
    ; Metrópolis criterium.
    ;-----------------------------------------
    fnew=call_function(funcname,xnew)
    if fnew le f then begin
      
      flag_aceptancia=1

    endif else begin

      ; This is used to prevent that the exponential cause overflow
      ; with (fnew-f)/Ta ~ +20. Actually, it can be used until 50 without
      ; problems but it doesn't make sense because 1./(1+exp(+20)) ~ 0.
      test=(fnew-f)/Ta	; como fnew > f => test > 0 sempre.
      Pa= test gt 20. ? 0d : 1./(1.+exp(test))
;       print,'test=',test,'  Pa=',Pa

      ; It is accepted con prob Punif.
      Punif=randomu(seed2,/double)	;entre 0 e 1.
      if Punif le Pa then flag_aceptancia=1

    endelse
    ;-----------------------------------------
    
    ; If there is an acceptance save the data in a file.
    if flag_aceptancia eq 1 then begin

    ; Se é o mellor f garda os parámetros.
      if fnew lt fbest then begin
	fbest=fnew
	xbest=xnew
	nrean=0
      endif
      
      nacep=nacep+1
      ka=ka+1
      x=xnew
      f=fnew
;       print,'ka='+strtrim(ka,2)+' Ta='+strtrim(Ta,2)+ $
;       	' fnew='+strtrim(fnew,2)+' fbest='+strtrim(fbest,2)

      if keyword_set(evolfile) then begin
	printf,uni, $
	  strtrim(Ta,2)+' '+ $
	  strtrim(f,2)+' '+ $
	  strtrim(x[0],2)+' '+ $
	  string(x[1],f='(d14.6)')+' '+ $
	  strtrim(x[2],2)+' '+ $
	  strtrim(x[3],2)+' '+ $
	  strtrim(x[4],2)+' '+ $
	  strtrim(x[5],2)+' '+ $
	  strtrim(x[6],2)
      endif
      
    endif

    ; Following Na acceptances do a reannealing.
    if nacep ge Na then begin

      print,'Re-annealing...'

      ; Compute the sensibilities s.
      s=replicate(0.0d,npar)
      for g=0,npar-1 do begin

	; Compute only for the fitted parameters.
	if fitparam[g] eq 1 then begin
	  
	  ee=replicate(0.,npar)
	  ee[g]=delta[g]
	  fbestdelta=call_function(funcname,xbest+ee)
	  s[g]=abs((fbestdelta-fbest)/delta[g])
	  
	endif
      endfor
      
      ; This is to avoid s=0 in denominator.
      ind0=where(s eq 0.0, count,complement=indno0)
      if count gt 0 then s[ind0]=min(s[indno0])
      
      smax=max(s[indfitted])
      
      ; Change the generating temperature and set kgen.
      Tgen[indfitted]=Tgen[indfitted]*(smax/s[indfitted])
      kgen[indfitted]=(alog(Tgen0[indfitted]/Tgen[indfitted])/c)^double(nfitted)
      kgen[indfitted]=abs(kgen[indfitted])

      ; Change the acceptance temperature and set ka.
      Ta0=f
      Ta=fbest
      ka=(alog(Ta0/Ta)/c)^double(nfitted)
      
      ; --------------------CHECK
;       print,'smax/s=',smax/s
;       print,'kgen=',kgen
      ;print,'ka=',ka,'  kgen=',kgen
      kk=where(finite(Tgen[indfitted]) eq 0,count)
      if count gt 0 then stop
      ; --------------------CHECK

      ; Reset counters.
      nacep=0UL
      nrean=nrean+1
      ;print,'nrean=',nrean
      
    endif
  endfor	; End of loop of generated points.

  ; Print the best values found for this acceptance temperature.
  print,"Ta="+string(Ta,format='(E10.4)')+ $
    " param=["+strtrim(xbest[0],2)+"," $
    +string(xbest[1],format='(d12.4)')+"," $
    +string(xbest[2],format='(d8.6)')+"," $
    +string(xbest[3],format='(d7.3)')+"," $
    +string(xbest[4],format='(d7.3)')+"," $
    +string(xbest[5],format='(d7.3)')+"," $
    +string(xbest[6],format='(d7.3)') $
    +"] chi^2="+strtrim(fbest,2)
  
  ; Reduction of generation temperature.
  kgen[indfitted]=kgen[indfitted]+1
  Tgen[indfitted]=Tgen0[indfitted]*exp(-c*kgen[indfitted]^(1d /nfitted))
  
  ; Reduction of acceptance temperature.
  ka=ka+1
  Ta=Ta0*exp(-c*ka^(1./nfitted))
  
  ; Place fbest at the end of fbestlist.
  fbestlist=[fbestlist[1:Neps-1],fbest]
  
  ; Check that the last Neps values of fbestlist are less than eps.
  diff=abs(fbestlist-shift(fbestlist,1))
  ind=where(diff lt eps,count)
  if count eq Neps then begin
    if fbest lt chislimit then break $	; Termination.
    else Ta=Ta0
  endif

  ; Ends if the number of reannealings with no improvemenent in fbest is equal to Nterm.
  if nrean ge Nterm then begin
    print,'Maximum number of reannealings reached.'
    break
  endif

endrep until 0	; End of the annealing loop.

; Close de evolution file if it was opened.
if keyword_set(evolfile) then begin
  close,uni
  FREE_LUN,uni
  print,'Parameters evolution file: '+evolfile
endif

; ------------------Compute the uncertainties-----------------------

; This is done by computing the Fisher matrix in the best point found.
; Usually this method underestimates the uncertainties, but it is only for
; a preliminary guess. If a detailed computation is needed another method
; (such as the MCMC) must be applied. Also, this guess can be used
; to run the MCMC near the solution.

; Fisher matrix.
FF=dblarr(nfitted,nfitted)

n1=0	; Counter
for g1=0,npar-1 do begin

  if fitparam[g1] eq 1 then begin
  
    ee1=replicate(0.,npar)
    ee1[g1]=delta[g1]
    
    n2=0
    for g2=0,npar-1 do begin

      if fitparam[g2] eq 1 then begin
        
	if g1 eq g2 then begin

	  fm=call_function(funcname,xbest-ee1)
	  fp=call_function(funcname,xbest+ee1)
	  ddf=(fp-2.*fbest+fm)/delta[g1]^2.

	endif else begin

	  ee2=replicate(0.,npar)
	  ee2[g2]=delta[g2]
	  fpp=call_function(funcname,xbest+ee1+ee2)
	  fpm=call_function(funcname,xbest+ee1-ee2)
	  fmp=call_function(funcname,xbest-ee1+ee2)
	  fmm=call_function(funcname,xbest-ee1-ee2)
	  ddf=(fpp-fpm-fmp+fmm)/(4.*delta[g1]*delta[g2])

	endelse
; 	print,'[g1, g2]=[',g1,' ,',g2,'], ddf=',ddf
	FF[n1,n2]=0.5*ddf
	
	n2=n2+1
      endif
    endfor
    
    n1=n1+1
  endif
endfor
; print,'FF='
; print,FF

;detF=determ(FF,/double)

; Covariance matrix.
cov=invert(FF,/double)
; print,'cov='
; print,cov

; Computation of the variances (elements in the diagonal).
; The abs inside the sqrt() is to avoid negative variances.
diag=indgen(nfitted)
sxbest=replicate(0.,npar)
sxbest[indfitted]=sqrt(abs(cov[diag,diag]))
; print,'cov[i,i]='
; print,cov[diag,diag]
; print,'sigmas='
; print,sxbest

; Time vector for the fitted curve.
tini=min(tt)
tfin=max(tt)

; Save the parameters and their uncertainties in a file.
writeparamfile,outfile,npar,xbest,sxbest,replicate(0.0d,npar),fbest,nobs1,nobs2,tfin-tini

; This part was extracted to an external program to re-run it separately if needed.
if keyword_set(physics) then $
  computeRVparams,readparfile=outfile,rvdatafile1=rvfile1,rvdatafile2=rvfile2,/fase

PRINT,'Processed in '+strtrim(SYSTIME(1)-tprocini,2)+' seconds.'

print,'FIN.'
end
