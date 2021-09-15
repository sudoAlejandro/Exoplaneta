;+
; NAME:
;	computeRVparams
;
; PURPOSE:
;	Compute the physical parameters for a spectroscopic binary.
;
; EXPLANATION:
;	This procedure can do several things: it computes the parameters and their
;	uncertainties from the keplerian parameters fitted, it plots the RV curve
;	with the residuals, and it can build a latex table with the parameters fitted.
;	It need the RVlib.pro file with the procedures to compute the radial velocities.
;	Also it need the pxperfect.pro file to produce pretty ps files.
;	
; CALLING SEQUENCE:
;	computeRVparams,readparfile=file,rvdatafile1=rvdatafile1,rvdatafile2=rvdatafile2
;		[,latextable=latextable][,\fase][,\ps]
;
; INPUTS:
;	readparfile = name of the parameters file as a result of rvfit (string).
;	rvdatafile1 = name of the file containing the primary radial velocities (string).
;	rvdatafile2 = name of the file containing the secondary radial velocities (string).
;	latextable = if set, name of the output .tex file with the results table (string).
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;	/fase = plot the RV data in phase with the period and Tp.
;	/ps = save a postscript file with the RV plot.
;
; OUTPUTS:
;	Values for physical parameters of the system, for a single-line or double-line
;	binary.
;	Plot with the RV curve and the residuals for the computed curve given by
;	the parameters in readparfile.
;	If set, .tex file with the parameters for publication.
;
; EXAMPLES:
;	Compute the physical parameters for the eclipsing binary LV Her and plot
;	the RV curve in phase with the period:
; 	computeRVparams,readparfile='LVHer.out',rvdatafile1='LVHer_RV1.dat',rvdatafile2='LVHer_RV2.dat',/fase
;
;	Compute the physical parameters for the exoplanet HD 37605, plot the RV curve
;	folded in phase and save the latex table with the results:
; 	computeRVparams,readparfile='HD37605.out',rvdatafile1='HD37605_RV.dat',/fase,/latextable
;
; REFERENCES:
;	
; MODIFICATION HISTORY:
;       Escrito por Ramon Iglesias Marzoa, Dep. Astrofisica (ULL), 2013-03-09.
;	Última modificación: 2015-02-20.


@RVlib

pro computeRVparams, $
  readparfile=readparfile, $
  rvdatafile1=rvdatafile1, $
  rvdatafile2=rvdatafile2, $
  fase=fase, $
  latextable=latextable, $
  ps=ps

; Number of parameters.
nx=7

; Read the parameters.
usage='computeRVparams,readparfile=file,rvdatafile1=file,rvdatafile2=file[,/fase][,/latextable]'
if ~keyword_set(readparfile) then begin
  print,'You must set a valid parameters file:'
  print,usage
  retall
endif else begin

  readparamfile,readparfile,nx,params,sinf,ssup,chis,nobs1,nobs2,dt
  P=params[0]
  Tp=params[1]
  e=params[2]
  omega=params[3]
  gamma=params[4]
  K1=params[5]
  K2=params[6]
  
  ; Check wether the uncertainties are symmetric.
  sel=where(ssup eq 0., count)
  if count eq nx then begin
    
    ; If all the upper uncert. are null then the actual uncert. are the lower.
    sP=sinf[0]
    sTp=sinf[1]
    se=sinf[2]
    somega=sinf[3]
    sgamma=sinf[4]
    sK1=sinf[5]
    sK2=sinf[6]

  endif else begin
    
    ; If the upper uncert. are not null then simmetrize them.
    sP=(sinf[0]+ssup[0])/2.
    sTp=(sinf[1]+ssup[1])/2.
    se=(sinf[2]+ssup[2])/2.
    somega=(sinf[3]+ssup[3])/2.
    sgamma=(sinf[4]+ssup[4])/2.
    sK1=(sinf[5]+ssup[5])/2.
    sK2=(sinf[6]+ssup[6])/2.
    
  endelse
endelse

tbase=2450000.0d

; Simbolo = círculo.
A=findgen(17*(!pi*2/16.))
A=[A,A[0]]
usersym,cos(A),sin(A),/fill
tsimbol=1.2	; tamaño

; This is for the plots.
if keyword_set(ps) then begin
  ; Save original graphic variables and create new ones.
  p_old=!p
  x_old=!x
  y_old=!y

  !p.thick=2
  !p.charthick=2
  !x.thick=2
  !y.thick=2
  
  ; Open postscript file
  psfile=readparfile+'.ps'
  ;read,prompt='Type postcript file name: ',psfile
  dummy=pxperfect()
  set_plot,'ps'
  device,_extra=dummy,filename=psfile,/color
endif else begin
;   wysize=500
;   wxsize=wysize*1.618	; golden rectangle.
  wxsize=640
  wysize=700
  device,retain=2
  window,0,xsize=wxsize,ysize=wysize
endelse
loadct,13
!p.charsize=1.7


if nobs2 gt 0 then begin

  ; This is for double-line binaries.

  if ~keyword_set(rvdatafile1) or ~keyword_set(rvdatafile2) then begin
    print,'You need to use valid data files.'
    retall
  endif
 
  ; Data reading.
  readcol,rvdatafile1,t1,rv1,srv1,f='d,f,f',comment='#',/silent
  readcol,rvdatafile2,t2,rv2,srv2,f='d,f,f',comment='#',/silent

  print
  print,'Result:'
  print,'------------------------------------------------'
  print,'              Fitted parameters'
  print,'------------------------------------------------'
  print,'P (d)              = '+strtrim(P,2)+' +/- '+strtrim(sP,2)
  print,'Tp (HJD/BJD)       = '+strtrim(string(Tp,f='(f13.5)'),2)+' +/- '+strtrim(sTP,2)
  print,'e                  = '+strtrim(e,2)+' +/- '+strtrim(se,2)
  print,'omega (deg)        = '+strtrim(omega,2)+' +/- '+strtrim(somega,2)
  print,'gamma (km/s)       = '+strtrim(gamma,2)+' +/- '+strtrim(sgamma,2)
  print,'K1 (km/s)          = '+strtrim(K1,2)+' +/- '+strtrim(sK1,2)
  print,'K2 (km/s)          = '+strtrim(K2,2)+' +/- '+strtrim(sK2,2)

  M1sini3=1.0361e-7*(1.-e^2)^(1.5)*(K1+K2)^2*K2*P	; Msun
  M2sini3=1.0361e-7*(1.-e^2)^(1.5)*(K1+K2)^2*K1*P	; Msun
  a1sini=13751.*sqrt(1.-e^2)*K1*P	; km
  a2sini=13751.*sqrt(1.-e^2)*K2*P	; km
  asini=a1sini+a2sini				; km
  q=M2sini3/M1sini3

  sM1sini3=1.0361e-7*sqrt( (3*(K1+K2)^2*K2*P*e)^2*(1-e^2)*se^2 + $
    (2*P*(K1+K2)*K2)^2*(1-e^2)^3*sK1^2 + $
    (2*(K1+K2)*K2+(K1+K2)^2)^2*P^2*(1-e^2)^3*sK2^2 + $
    (K1+K2)^4*K2^2*(1-e^2)^3*sP^2)
  sM2sini3=1.0361e-7*sqrt( (3*(K1+K2)^2*K1*P*e)^2*(1-e^2)*se^2 + $
    (2*(K1+K2)*K1+(K1+K2)^2)^2*P^2*(1-e^2)^3*sK1^2 + $
    (2*P*(K1+K2)*K1)^2*(1-e^2)^3*sK2^2 + $
    (K1+K2)^4*K1^2*(1-e^2)^3*sP^2)
  sa1sini=13751.*sqrt((1-e^2)*( (K1*P*e/(1-e^2))^2*se^2 + P^2*sK1^2 + K1^2*sP^2))
  sa2sini=13751.*sqrt((1-e^2)*( (K2*P*e/(1-e^2))^2*se^2 + P^2*sK2^2 + K1^2*sP^2))
  sasini=sqrt(sa1sini^2+sa2sini^2)
  sq=sqrt((sK1/K2)^2+(K1*sK2/K2^2)^2)

  print,'------------------------------------------------'
  print,'              Derived quantities'
  print,'------------------------------------------------'
  print,'M1sin(i)^3 (Msun)  = '+strtrim(M1sini3,2)+' +/- '+strtrim(sM1sini3,2)
  print,'M2sin(i)^3 (Msun)  = '+strtrim(M2sini3,2)+' +/- '+strtrim(sM2sini3,2)
  print,'q = M2/M1          = '+strtrim(q,2)+' +/- '+strtrim(sq,2)
  print,'a1sin(i) (10^6 km) = '+strtrim(a1sini/1e6,2)+' +/- '+strtrim(sa1sini/1e6,2)
  print,'         (Rsun)    = '+strtrim(a1sini*0.019758/13751.,2)+' +/- '+strtrim(sa1sini*0.019758/13751.,2)
  print,'a2sin(i) (10^6 km) = '+strtrim(a2sini/1e6,2)+' +/- '+strtrim(sa2sini/1e6,2)
  print,'         (Rsun)    = '+strtrim(a2sini*0.019758/13751.,2)+' +/- '+strtrim(sa2sini*0.019758/13751.,2)
  print,'asin(i) (10^6 km)  = '+strtrim(asini/1e6,2)+' +/- '+strtrim(sasini/1e6,2)
  print,'        (Rsun)     = '+strtrim(asini*0.019758/13751.,2)+' +/- '+strtrim(sasini*0.019758/13751.,2)

  ; Time array for the fitted curve.
  tt=[t1,t2]
  rrvv=[rv1,rv2]
  tini=min(tt)
  tfin=max(tt)
  np=(tfin-tini)/P*100.	; 100 points in each observed period.
  tfit=tini+(tfin-tini)*findgen(np)/(np-1)
  ;print,tfit

  ; Array with the fitted radial velocities.
  vfit1=radvel2(tfit,[P,Tp,e,omega,gamma,K1])
  vfit2=radvel2(tfit,[P,Tp,e,omega+180.,gamma,K2])

  ; Compute the residuals.
  vinterpol1=interpol(vfit1,tfit,t1)
  OC1=rv1-vinterpol1
  vinterpol2=interpol(vfit2,tfit,t2)
  OC2=rv2-vinterpol2

  rms1=sqrt(total(OC1^2)/nobs1)
  rms2=sqrt(total(OC2^2)/nobs2)

  print,'------------------------------------------------'
  print,'              Other quantities'
  print,'------------------------------------------------'
  print,'chi^2              = '+strtrim(chis,2)
  print,'Nobs (primary)     = '+strtrim(nobs1,2)
  print,'Nobs (secondary)   = '+strtrim(nobs2,2)
  print,'Time span (days)   = '+strtrim(dt,2)
  print,'rms1 (km/s)        = '+strtrim(rms1,2)
  print,'rms2 (km/s)        = '+strtrim(rms2,2)

  ;-----------------LaTeX table--------------------

  if keyword_set(latextable) then begin
    latexfile=readparfile+'.tex'
    openw,uni,latexfile,/GET_LUN,width=160
    
    printf,uni,'\begin{center}'
    printf,uni,'{\scriptsize'
    printf,uni,'\begin{table}[h]'
    printf,uni,'\begin{tabular}{lr}'
    printf,uni,'\multicolumn{2}{c}{TITLE OF THE TABLE}\\'
    printf,uni,'\hline'
    printf,uni,'\hline'
    printf,uni,'Parameter		&Value\\'
    printf,uni,'\hline'
    printf,uni,'\multicolumn{2}{c}{Adjusted Quantities}\\'
    printf,uni,'\hline'
    printf,uni,'$P$ (d)		&'+strtrim(P,2)+' $\pm$ '+strtrim(sP,2)+'\\'
    printf,uni,'$T_p$ (HJD)		&'+strtrim(string(Tp,f='(f13.5)'),2)+' $\pm$ '+strtrim(sTP,2)+'\\'
    printf,uni,'$e$			&'+strtrim(e,2)+' $\pm$ '+strtrim(se,2)+'\\'
    printf,uni,'$\omega$ (deg)		&'+strtrim(omega,2)+' $\pm$ '+strtrim(somega,2)+'\\'
    printf,uni,'$\gamma$ (km/s)	&'+strtrim(gamma,2)+' $\pm$ '+strtrim(sgamma,2)+'\\'
    printf,uni,'$K_1$ (km/s)		&'+strtrim(K1,2)+' $\pm$ '+strtrim(sK1,2)+'\\'
    printf,uni,'$K_2$ (km/s)		&'+strtrim(K2,2)+' $\pm$ '+strtrim(sK2,2)+'\\'
    printf,uni,'\hline'
    printf,uni,'\multicolumn{2}{c}{Derived Quantities}\\'
    printf,uni,'\hline'
    printf,uni,'$M_1\sin ^3i$ ($M_\odot$)	&'+strtrim(M1sini3,2)+' $\pm$ '+strtrim(sM1sini3,2)+'\\'
    printf,uni,'$M_2\sin ^3i$ ($M_\odot$)	&'+strtrim(M2sini3,2)+' $\pm$ '+strtrim(sM2sini3,2)+'\\'
    printf,uni,'$q = M_2/M_1$		&'+strtrim(q,2)+' $\pm$ '+strtrim(sq,2)+'\\'
    printf,uni,'$a_1\sin i$ ($10^6$ km)	&'+strtrim(a1sini/1e6,2)+' $\pm$ '+strtrim(sa1sini/1e6,2)+'\\'
    printf,uni,'$a_2\sin i$ ($10^6$ km)	&'+strtrim(a2sini/1e6,2)+' $\pm$ '+strtrim(sa2sini/1e6,2)+'\\'
    printf,uni,'$a  \sin i$ ($10^6$ km)	&'+strtrim(asini/1e6,2)+' $\pm$ '+strtrim(sasini/1e6,2)+'\\'
    printf,uni,'\hline'
    printf,uni,'\multicolumn{2}{c}{Other Quantities}\\'
    printf,uni,'\hline'
    printf,uni,'$\chi^2$			&'+strtrim(chis,2)+'\\'
    printf,uni,'$N_{obs}$ (primary)		&'+strtrim(nobs1,2)+'\\'
    printf,uni,'$N_{obs}$ (secondary)	&'+strtrim(nobs2,2)+'\\'
    printf,uni,'Time span (days)		&'+strtrim(dt,2)+'\\'
    printf,uni,'$rms_1$ (km/s)		&'+strtrim(rms1,2)+'\\'
    printf,uni,'$rms_2$ (km/s)		&'+strtrim(rms2,2)+'\\'
    printf,uni,'\hline'
    printf,uni,'\end{tabular}'
    printf,uni,'\caption{\footnotesize $^a$ Parameter fixed beforehand.}'
    printf,uni,'\label{table:test2}'
    printf,uni,'\end{table}'
    printf,uni,'}'
    printf,uni,'\end{center}'
    
    close,uni
    FREE_LUN,uni
    print,'LaTeX file: '+latexfile
  
  endif
  
  ;-----------------Plots--------------------

  rangot=[tini,tfin]-tbase
  rangorv=[min([rrvv,vfit1,vfit2]),max([rrvv,vfit1,vfit2])]
  rangoOC=[min([OC1,OC2]),max([OC1,OC2])]
  
  pos1=[0.15,	0.535,			0.96,	0.95]
  pos2=[0.15,	0.11+0.2025+0.01,	0.96,	0.535-0.01]
  pos3=[0.15,	0.11, 			0.96,	0.11+0.2025]

  if keyword_set(fase) then begin

    ; Fases.
    ciclo=(t1-Tp)/P	; ciclo e o instante da observacion medido en unidades de periodo
    fase1=ciclo-floor(ciclo)		; floor quedase coa parte enteira do ciclo
    ordenobs1=sort(fase1)
    fase1=fase1(ordenobs1)

    ciclo=(t2-Tp)/P	; ciclo e o instante da observacion medido en unidades de periodo
    fase2=ciclo-floor(ciclo)		; floor quedase coa parte enteira do ciclo
    ordenobs2=sort(fase2)
    fase2=fase2(ordenobs2)

    ciclo=(tfit-Tp)/P	; ciclo e o instante da observacion medido en unidades de periodo
    fase=ciclo-floor(ciclo)		; floor quedase coa parte enteira do ciclo
    ordenfit=sort(fase)
    fasefit=fase(ordenfit)

    wtitle='RV fit and residuals'
    xx1=fase1
    xx2=fase2
    yy1=rv1[ordenobs1]
    yy2=rv2[ordenobs2]
    rr1=OC1[ordenobs1]
    rr2=OC2[ordenobs2]
    syy1=srv1[ordenobs1]
    syy2=srv2[ordenobs2]
    xxfit=fasefit
    yyfit1=vfit1[ordenfit]
    yyfit2=vfit2[ordenfit]
    xxrange=[0.,1.]
    yyrange=[min(rrvv),max(rrvv)]
    yyrange=rangorv
;     rrrange=[min([OC1,OC2]),max([OC1,OC2])]
    rrrange=max(abs([OC1,OC2]))*[-1.,1.]
    xxtitle='Phase'
    yytitle='RV (km/s)'
    rrtitle1='(O-C)!D1!N (km/s)'
    rrtitle2='(O-C)!D2!N (km/s)'
    
  endif else begin

    wtitle='RV fit and residuals'
    xx1=t1-tbase
    xx2=t2-tbase
    yy1=rv1
    yy2=rv2
    rr1=OC1
    rr2=OC2
    syy1=srv1
    syy2=srv2
    xxfit=tfit-tbase
    yyfit1=vfit1
    yyfit2=vfit2
    xxrange=[tini,tfin]-tbase
;     yyrange=[min(rrvv),max(rrvv)]
    yyrange=rangorv
;     rrrange=[min([OC1,OC2]),max([OC1,OC2])]
    rrrange=max(abs([OC1,OC2]))*[-1.,1.]
    xxtitle='HJD-'+strtrim(tbase,2)
    yytitle='RV (km/s)'
    rrtitle1='(O-C)!D1!N (km/s)'
    rrtitle2='(O-C)!D2!N (km/s)'
    
  endelse

  ; To plot curves with only 1 RV observed point in one component
  ; (usually the secondary).
  if n_elements(xx1) eq 1 then begin
    xx1=replicate(xx1,2)
    yy1=replicate(yy1,2)
    syy1=replicate(syy1,2)
    rr1=replicate(rr1,2)
  endif
  if n_elements(xx2) eq 1 then begin
    xx2=replicate(xx2,2)
    yy2=replicate(yy2,2)
    syy2=replicate(syy2,2)
    rr2=replicate(rr2,2)
  endif
  
  ; Plots with the measurements and the fit.
  ploterror,xx1,yy1,syy1,psym=8,/nohat, $
    position=pos1, $
    xrange=xxrange,XTICKFORMAT="(A1)", $
    ytitle=yytitle,yrange=yyrange
  oplot,xxfit,yyfit1,psym=0,color=2000

  usersym,cos(A),sin(A)
  oploterror,xx2,yy2,syy2,psym=8,/nohat
  oplot,xxfit,yyfit2,psym=0,linestyle=2,color=2000

  plots,!x.crange,gamma*[1.,1.],color=2000,linestyle=1	; gamma line.

  ; Plots with residuals.
  usersym,cos(A),sin(A),/fill
  ploterror,xx1,rr1,syy1,psym=8,/nohat,/noerase, $
    position=pos2, $
    xrange=xrange,XTICKFORMAT="(A1)", $
    ytitle=rrtitle1,yrange=rrrange,ystyle=1,ytickv=rr1[0]
  plots,!x.crange,[0.,0.],color=2000

  usersym,cos(A),sin(A)
  ploterror,xx2,rr2,syy2,psym=8,/nohat,/noerase, $
    position=pos3, $
    xtitle=xxtitle,xrange=xxrange, $
    ytitle=rrtitle2,yrange=rrrange,ystyle=1,ytickv=rr2[0]
  plots,!x.crange,[0.,0.],color=2000,linestyle=2
  
endif else begin

  ; This is for single-line binaries.

  if ~keyword_set(rvdatafile1) then begin
    print,'You need to use valid data files.'
    retall
  endif
 
  ; Data reading.
  readcol,rvdatafile1,t1,rv1,srv1,f='d,f,f',comment='#',/silent

  print
  print,'Result:'
  print,'------------------------------------------------'
  print,'              Fitted parameters'
  print,'------------------------------------------------'
  print,'P (d)              = '+strtrim(P,2)+' +/- '+strtrim(sP,2)
  print,'Tp (HJD/BJD)       = '+strtrim(string(Tp,f='(f13.5)'),2)+' +/- '+strtrim(sTP,2)
  print,'e                  = '+strtrim(e,2)+' +/- '+strtrim(se,2)
  print,'omega (deg)        = '+strtrim(omega,2)+' +/- '+strtrim(somega,2)
  print,'gamma (km/s)       = '+strtrim(gamma,2)+' +/- '+strtrim(sgamma,2)
  print,'K1 (km/s)          = '+strtrim(K1,2)+' +/- '+strtrim(sK1,2)
  
  a1sini=13751.*sqrt(1.-e^2)*K1*P	; Km 
  f=1.0361e-7*(1-e^2)^(3./2.)*K1^3*P

  sa1sini=13751.*sqrt((1-e^2)*( (K1*P*e/(1-e^2))^2*se^2 + P^2*sK1^2 + K1^2*sP^2))
  sf=1.0361e-7*sqrt((9*K1^6*P^2*(1-e^2))*se^2 + (9*(1-e^2)^3*K1^4*P^2)*sK1^2 + ((1-e^2)^3*K1^6)*sP^2)

  print,'------------------------------------------------'
  print,'              Derived quantities'
  print,'------------------------------------------------'
  print,'a1sin(i) (10^6 km) = '+strtrim(a1sini/1e6,2)+' +/- '+strtrim(sa1sini/1e6,2)
  print,'         (Rsun)    = '+strtrim(a1sini*0.019758/13751.,2)+' +/- '+strtrim(sa1sini*0.019758/13751.,2)
  print,'f(m1,m2) (Msun)    = '+strtrim(f,2)+' +/- '+strtrim(sf,2)

  ; Time array for the fitted curve.
  tini=min(t1)
  tfin=max(t1)
  np=(tfin-tini)/P*100.	; 100 points in each observed period.
  tfit=tini+(tfin-tini)*findgen(np)/(np-1)
  ;print,tfit

  ; Array with the fitted radial velocities.
;   if e lt 1e-4 then begin
;     e=0.
;     omega=0.
;   endif
  vfit1=radvel2(tfit,[P,Tp,e,omega,gamma,K1])

  ; Compute the residuals.
  vinterpol1=interpol(vfit1,tfit,t1)
  OC1=rv1-vinterpol1
  rms1=sqrt(total(OC1^2)/nobs1)

  print,'------------------------------------------------'
  print,'              Other quantities'
  print,'------------------------------------------------'
  print,'chi^2              = '+strtrim(chis,2)
  print,'Nobs (primary)     = '+strtrim(nobs1,2)
  print,'Time span (days)   = '+strtrim(dt,2)
  print,'rms1 (km/s)        = '+strtrim(rms1,2)
  
  ;-----------------LaTeX table--------------------

  if keyword_set(latextable) then begin
  
    latexfile=readparfile+'.tex'
    openw,uni,latexfile,/GET_LUN,width=160
    
    printf,uni,'\begin{center}'
    printf,uni,'{\scriptsize'
    printf,uni,'\begin{table}[h]'
    printf,uni,'\begin{tabular}{lr}'
    printf,uni,'\multicolumn{2}{c}{TITLE OF THE TABLE}\\'
    printf,uni,'\hline'
    printf,uni,'\hline'
    printf,uni,'Parameter			&Value\\'
    printf,uni,'\hline'
    printf,uni,'\multicolumn{2}{c}{Adjusted Quantities}\\'
    printf,uni,'\hline'
    printf,uni,'$P$ (d)		&'+strtrim(P,2)+' $\pm$ '+strtrim(sP,2)+'\\'
    printf,uni,'$T_p$ (HJD)		&'+strtrim(string(Tp,f='(f13.5)'),2)+' $\pm$ '+strtrim(sTP,2)+'\\'
    printf,uni,'$e$			&'+strtrim(e,2)+' $\pm$ '+strtrim(se,2)+'\\'
    printf,uni,'$\omega$ (deg)		&'+strtrim(omega,2)+' $\pm$ '+strtrim(somega,2)+'\\'
    printf,uni,'$\gamma$ (km/s)	&'+strtrim(gamma,2)+' $\pm$ '+strtrim(sgamma,2)+'\\'
    printf,uni,'$K_1$ (km/s)		&'+strtrim(K1,2)+' $\pm$ '+strtrim(sK1,2)+'\\'
    printf,uni,'\hline'
    printf,uni,'\multicolumn{2}{c}{Derived Quantities}\\'
    printf,uni,'\hline'
    printf,uni,'$a_1\sin i$ ($10^6$ km)	&'+strtrim(a1sini/1e6,2)+' $\pm$ '+strtrim(sa1sini/1e6,2)+'\\'
    printf,uni,'$f(m_1,m_2)$ ($M_\odot$)	&'+strtrim(f,2)+' $\pm$ '+strtrim(sf,2)+'\\'
    printf,uni,'\hline'
    printf,uni,'\multicolumn{2}{c}{Other Quantities}\\'
    printf,uni,'\hline'
    printf,uni,'$\chi^2$		&'+strtrim(chis,2)+'\\'
    printf,uni,'$N_{obs}$ (primary)	&'+strtrim(nobs1,2)+'\\'
    printf,uni,'Time span (days)	&'+strtrim(dt,2)+'\\'
    printf,uni,'$rms_1$ (km/s)	&'+strtrim(rms1,2)+'\\'
    printf,uni,'\hline'
    printf,uni,'\end{tabular}'
    printf,uni,'\caption{\footnotesize $^a$ Parameter fixed beforehand.}'
    printf,uni,'\label{table:test4}'
    printf,uni,'\end{table}'
    printf,uni,'}'
    printf,uni,'\end{center}'
    
    close,uni
    FREE_LUN,uni
    print,'LaTeX file: '+latexfile
  
  endif
  
  ;-----------------Plots--------------------

  pos1=[0.15,	0.11+0.2025+0.01,	0.96,	0.95]
  pos2=[0.15,	0.11,			0.96,	0.11+0.2025]

  if keyword_set(fase) then begin

    ; Fases.
    ciclo=(t1-Tp)/P	; ciclo e o instante da observacion medido en unidades de periodo
    fase1=ciclo-floor(ciclo)		; floor quedase coa parte enteira do ciclo
    ordenobs1=sort(fase1)
    fase1=fase1(ordenobs1)

    ciclo=(tfit-Tp)/P	; ciclo e o instante da observacion medido en unidades de periodo
    fasefit=ciclo-floor(ciclo)		; floor quedase coa parte enteira do ciclo
    ordenfit=sort(fasefit)
    fasefit=fasefit(ordenfit)

    wtitle='RV fit and residuals (phased)'
    xx1=fase1
    yy1=rv1[ordenobs1]
    rr1=OC1[ordenobs1]
    syy1=srv1[ordenobs1]
    xxfit=fasefit
    yyfit=vfit1[ordenfit]
    xrange=[0.,1.]
    yrange=[min([rv1,vfit1]),max([rv1,vfit1])]
    rrange=max(abs(OC1))*[-1.,1.]
    xtitle='Phase'
    ytitle1='RV (km/s)'
    rrtitle='(O-C)!D!N (km/s)'

  endif else begin
    
    wtitle='RV fit and residuals'
    xx1=t1-tbase
    yy1=rv1
    rr1=OC1
    syy1=srv1
    xxfit=tfit-tbase
    yyfit=vfit1
    xrange=[min(t1),max(t1)]-tbase
    yrange=[min([rv1,vfit1]),max([rv1,rvfit1])]
    rrange=max(abs(OC1))*[-1.,1.]
    xtitle='HJD-'+strtrim(tbase,2)
    ytitle1='RV (km/s)'
    rrtitle='(O-C)!D!N (km/s)'
 
  endelse

  ; Plots with the measurements and the fit.
  ploterror,xx1,yy1,syy1,psym=8,/nohat, $
    position=pos1, $
    xrange=xrange,XTICKFORMAT="(A1)", $
    ytitle=ytitle1,ystyle=1,ytickv=rr1[0],yrange=yrange
  oplot,xxfit,yyfit,psym=0,color=2000

  plots,!x.crange,gamma*[1.,1.],color=2000,linestyle=1	; gamma line.

  ; Plots with residuals.
  ploterror,xx1,rr1,syy1,psym=8,/nohat,/noerase, $
    position=pos2, $
    xtitle=xtitle,xrange=xrange, $
    ytitle=rrtitle,yrange=rrange,ystyle=1,ytickv=rr1[0]

  plots,!x.crange,[0.,0.],color=2000
  
endelse

; Close postcript file and restore graphic variables.
if keyword_set(ps) then begin
  device,/close
  print,'Output: '+psfile
  set_plot,'x'
  !p=p_old
  !x=x_old
  !y=y_old
endif

end
