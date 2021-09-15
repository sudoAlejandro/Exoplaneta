;print,'baseado en Navarro&Perfors con xeneración q(x) gausiana.'
; http://health.adelaide.edu.au/psychology/ccs/docs/ccs-class/technote_metropolishastings.pdf

; This program asumes that you have run rvfit and you got a parameters file from it.
;MCMCerrors,'LVHer_param.out','LVHer_RV1.dat','LVHer_RV2.dat',outfile='LVHer_MCMC.out',nbins=30,C=0.9545,nsample=ulong(1e6)
;MCMCerrors,'LVHer_param.out','LVHer_RV1.dat',''

; RIM, 2013-06-14.
; Última modificación: 2014-11-04.

; Engade as funcións necesarias.
@RVlib

; Función de xeneración.
function gen,x,U,L,fitparam,sigma,gaus=gaus
nx=n_elements(x)
xnew=x

; Percorre as direccións xenerando un novo valor.
for g=0,nx-1 do begin

  ; Se este parámetro está fixado pasa ó seguinte.
  if fitparam[g] eq 1 then begin

    ; Xenera novos valores da función dentro do dominio.
    repeat begin
    
      if keyword_set(gaus) then begin
      
	; Distribución gausiana.
	rand=randomn(seed1,/double)	; rand distribuido de forma gausiana en torno a 0 con sigma=1.
	xnew[g]=x[g]+sigma[g]*rand	; Fai o cambio para centrar en xnew con anchura sigma en vez de 1.
      
      endif else begin
      
	; Distribución uniforme.
	rand=-1d0+2d0*randomu(seed1,/double)	; r entre -1 e 1.
	xnew[g]=x[g]+rand*(U[g]-L[g])    
      
      endelse
    endrep until xnew[g] ge L[g] and xnew[g] le U[g]
;     print,'g=',g,' rand=',rand, $
;       ' U[g]=',U[g],' L[g]=',L[g],' xnew[g]=',xnew[g]
  endif
endfor	; Fin do bucle de intentos en cada dirección.

return,xnew
end


pro MCMCerrors,bestparfile,rvdatafile1,rvdatafile2,outfile=outfile,nbins=nbins,C=C,nsample=nsample

common datos,t1,rv1,srv1,t2,rv2,srv2

; junk=check_math(/print)
; !except=2

usage='MCMCerrors,bestparfile,rvdatafile1,rvdatafile2,[,outfile=outfile][,nbins=nbins][,C=C][,nsample=nsample]'
if n_params() lt 3 then begin
  print,'Usage:'
  print,usageconfigfile
  retall
endif

; Set default parameters to run the program.
if ~keyword_set(outfile) then outfile='MCMCerrors.out'
if ~keyword_set(nbins) then nbins=50
if ~keyword_set(C) then C=0.6827
if ~keyword_set(nsample) then nsample=100000UL
; C=0.6827 for 68.3% confidence interval.
; C=0.9545 for 95.4% confidence interval.
; C=0.9973 for 99.7% confidence interval.

; Read the resulting keplerian parameters from ASA.
npar=7
readparamfile,bestparfile,npar,xbest,sl,sh,chis,nobs1,nobs2,dt

; If all the upper sigmas are null then the uncertainties are simmetrical and
; they are equal to the lower sigmas.
selection=where(sh eq 0.,count)
if count eq npar then sh=sl

; Rebuild fitparam using the computed uncertainties in bestparfile.
fitparam=replicate(1,npar)
indfit=where(sl eq 0.,count)
if count gt 0 then fitparam[indfit]=0

; Double-line flag to distinguish between double and single-line binaries.
if rvdatafile2 eq '' or nobs2 eq 0 then flagDL=0 else flagDL=1

; Lectura de datos.
readcol,rvdatafile1,t1,rv1,srv1,f='d,f,f',comment='#',/silent
;t1=t1+2400000d	; para os datos de LVHer

if flagDL eq 1 then begin
readcol,rvdatafile2,t2,rv2,srv2,f='d,f,f',comment='#',/silent
  funcname='chisqDL'
endif else begin
  fitparam[6]=0
  funcname='chisqSL'
endelse

; ------------------Dominio de cálculo------------------

; Build the computing domain.
nsig=10.
L=xbest-nsig*sl
U=xbest+nsig*sh

; Adjust the limits if they are outside the allowed values.
L[0]=max([L[0],1e-15])		; P -> [0.00...1, inf]
L[2]=max([L[2],0.0d])		; e -> [0, 0.999...]
U[2]=min([U[2],1.0d - 1.0e-9])	; e -> [0, 0.999...]
L[5]=max([L[5],0.])		; K1 -> [0, inf]
L[6]=max([L[6],0.])		; K2 -> [0, inf]

;------------------------------------------------
; Initialization of MCMC parameters.
;------------------------------------------------

x=xbest
nx=n_elements(x)
xnew=x
f=call_function(funcname,x)
burnin=0
lag=1

;------------------------------------------------
; Generate a suitable sigma for the gaussian of the generation function. It is
; considered good if the acceptance rate is ~25%.
;------------------------------------------------

acceptarget=0.25		; Target acceptance.
sigma=abs(U-L)/3.	; Initial test sigma.
nsamplesigma=1000	; Number or samples to measure the acceptance.

print,'Calibrating the sigma for target acceptance=',strtrim(acceptarget,2)
repeat begin
  naccep=0
  sigma=sigma/1.2
  
  x=xbest
  for j=0,nsamplesigma-1 do begin

    ; Xenera un novo valor.
    xnew=gen(x,U,L,fitparam,sigma,/gaus)
    fnew=call_function(funcname,xnew)
  
    uu=randomu(seed2,/double)	;entre 0 e 1.
    logr=0.5*(f-fnew)
    ; Para evitar underflow se logr << 0.
    if logr le -20. then logr=-20.
    ; Para evitar overflow, ver "Handbook of Markov Chain Monte Carlo", Brooks+ (ed), p 23.
    ; O segundo operando do or só se evalúa se o primeiro é 0, igual en IDL que en C.    
    if (logr ge 0.) or (uu lt exp(logr)) then begin
      x=xnew
      f=fnew
      naccep=naccep+1
    endif
  endfor

  acceprate=float(naccep)/float(nsamplesigma)
  print,strtrim(naccep,2),strtrim(acceprate,2),String(13b),format='("naccep=",A,"	acceprate=",A,A,$)'
endrep until acceprate gt acceptarget
print

; Reset the initial position to the parameters values.
x=xbest
f=call_function(funcname,xbest)

; Initial parameters.
nplot=nsample/50
jlast=0UL
markovchain=dblarr(nx,nsample)
naccep=0UL
nrun=0UL

; For plotting.
parname=['P','T!Dp!N','e','!7x!3','!7c!3','K!D1!N','K!D2!N']
wxsize=1500
wysize=700
window,0,xsize=wxsize,ysize=wysize,title='Parameter distribution'
!p.multi=[0,npar,2]

openw,uni,outfile,/GET_LUN,width=160

for j=0UL,nsample-1 do begin

  for k=0,lag-1 do begin
  
    ; Xenera un novo valor.
    xnew=gen(x,U,L,fitparam,sigma,/gaus)
    fnew=call_function(funcname,xnew)
;     printf,uni,1.,fnew,xnew

    uu=randomu(seed2,/double)	;entre 0 e 1.
    logr=0.5*(f-fnew)
    ; Para evitar underflow se logr << 0.
    if logr le -20. then logr=-20.
    ; Para evitar overflow, ver "Handbook of Markov Chain Monte Carlo", Brooks+ (ed), p 23.
    ; O segundo operando do "or" só se evalúa se o primeiro é 0, igual en IDL que en C.    
    if (logr ge 0.) or (uu lt exp(logr)) then begin
      x=xnew
      f=fnew
      accep=1
    endif else begin
      accep=0
    endelse
    nrun=nrun+1
    naccep=naccep+accep
    
  endfor ; Final del bucle de lag.
  markovchain[*,j]=x
  
  ; Save to file.
  printf,uni,string(x,f='(d12.7," ",d13.5," ",d10.8," ",d10.6," ",d12.6," ",d11.6," ",d11.6)')
  
  ; Plot the results.
  
  ; This is for plotting the last nplot chain data.
  if ((j mod nplot) eq 0UL and j gt 0) or (j eq nsample-1) then begin

    acceprate=float(naccep)/float(nrun)
    print,'nrun='+strtrim(nrun,2)
    print,'naccep='+strtrim(naccep,2)
    print,'acceprate='+strtrim(acceprate,2)
    print,'j=',strtrim(j,2)
    print,'----------------'
    
    ; Liña superior de histogramas.
    xrange=replicate(0.,2,nx)
    for g=0,nx-1 do begin
    
      if fitparam[g] eq 1 then begin

	h=histogram(markovchain[g,0:j],nbins=nbins,locations=locs)
	
	if g eq 1 then begin	; This is for Tp.
	  Tpbase=ulong(min(locs))
	  locs=locs-Tpbase
	  offset=Tpbase
; 	  title=parname[g]+'-'+string(Tpbase,format='(d13.5)')
	  title=parname[g]+'-'+strtrim(Tpbase,2)
	endif else begin
	  offset=0.
	  title=parname[g]
	endelse

	; Save the xrange for it use in the Markov chain plots in the bottom line.
	xrange[*,g]=[min(locs),max(locs)]
	
	; Plot histogram.
	plot,locs,h,psym=10,xrange=xrange[*,g],xstyle=1,charsize=2.0,title=title,xticks=2
	
	plots,xbest[g]*[1.,1.]-offset,!y.crange,color=2000	; raya vertical en el valor calculado.
	
	intervaloConf=ConfIntShort(h,locs,C,xbest[g])
	plots,intervaloConf[0]*[1.,1.],!y.crange,color=2000,linestyle=2
	;plots,intervaloConf[1]*[1.,1.],!y.crange,color=20000,linestyle=1
	plots,intervaloConf[2]*[1.,1.],!y.crange,color=2000,linestyle=2
	
	; This is for fitting a gaussian to the histogram.
;         hfit=gaussfit(locs,h,nterms=3,sigma=asigma)
;         oplot,locs,hfit,color=15020,linestyle=5,thick=3
      endif else begin
      
	plot,[0.],xrange=[0.,1.],charsize=2.0,title=parname[g],xticks=2,/nodata

      endelse
    endfor

    ; Bottom line with the Markov chain plots.
    for g=0,nx-1 do begin
      if fitparam[g] eq 1 then begin
      
	; This is for Tp.
	if g eq 1 then offset=Tpbase else offset=0.
	
	plot,markovchain[g,jlast:j]-offset,jlast+ulindgen(n_elements(markovchain[g,jlast:j])),charsize=2.0, $
	  xstyle=1,xrange=xrange[*,g],/ynozero,xticks=2
      endif else begin
	plot,[0.],xrange=[0.,1.],charsize=2.0,title=parname[g],xticks=2,/nodata
      endelse
    endfor
;--------------------------------------

    jlast=j
  endif

endfor
close,uni
FREE_LUN,uni
print,'Output file: '+outfile

!p.multi=0

end
