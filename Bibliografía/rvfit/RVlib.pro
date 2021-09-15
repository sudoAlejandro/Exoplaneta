; RIM, 2013-06-14.
; Last change: 2014-12-30.

;-------------------------------------------------------------------------------

; Write parameters file.
pro writeparamfile,filename,npar,param,low,high,chis,nobs1,nobs2,dt

openw,uni,filename,/GET_LUN
for g=0,npar-1 do printf,uni,param[g],low[g],high[g],format='(3(d-14.6,:," "))'
printf,uni,strtrim(chis,2)
printf,uni,strtrim(nobs1,2)
printf,uni,strtrim(nobs2,2)
printf,uni,strtrim(dt,2)
close,uni
FREE_LUN,uni
print,'Out file: '+filename

end
;-------------------------------------------------------------------------------

; Read parameters file.
pro readparamfile,filename,np,pars,lo,hi,chs,nn1,nn2,dt

; Read the parameters.
readcol,filename,all,f='d',comment='#',/silent
nlines=n_elements(all)
if nlines eq 11 then begin
  pars=double(all[0:6])
  chs=double(all[7])
  nn1=long(all[8])
  nn2=long(all[9])
  dt=double(all[10])

  ; Read the uncertainties.
  readcol,filename,lo,hi,f='x,d,d',comment='#',/silent
  
endif else begin
  print,'You must use a valid parameters file:'
  retall
endelse


end

;-------------------------------------------------------------------------------

; Solve the Kepler equation.
; M must be in radians.
function keplerec2,M,e
M=double(M)
e=double(e)

; Valor aproximado inicial para reducir ó mínimo as iteracións.
; Calcúlase usando un desenrolo en serie de MacLaurin de potencias
; de e usando E como parámetro, ver Apuntes de Astronomía, U.
; de Barcelona, ecuación 47.3, ou Heintz, p35.
;EE=M+e*sin(M)+e^2/2.*sin(2.*M)
EE=M+e*sin(M)+e^2/2.*sin(2.*M)+e^3/8.*(3*sin(3*M)-sin(M))

eps=1d-10	; Precission.
j=0	; Iteration counter.
repeat begin

  ; Newton method (Heintz, p 35).
  EE0=EE
  EE=EE0+(M+e*sin(EE0)-EE0)/(1.-e*cos(EE0))

  ; Kepler method.
  ; Ver Apuntes de Astronomía, U. de Barcelona, apartado 3.6.2
  ;EE0=EE
  ;EE=M+e*sin(EE)

  ;print,'EE=',EE,' j=',j
  j=j+1
endrep until max(abs(EE0-EE)) le eps
;print,'EE=',EE

return,EE
end

;-------------------------------------------------------------------------------

; Compute RV from the orbital parameters:
; P in days.
; Tp in Heliocentric Julian Days
; omega in degrees.
; K, gamma in anyone units (m/s or Km/s).
function radvel2,t,param

P=param[0]
Tp=param[1]
e=param[2]
omega=param[3]*!dpi/180.	; in rads.
gamma=param[4]
K=param[5]

; Mean anomaly, MOD 1. is to bring it into 0 and 2pi.
M=2.*!dpi*((t-Tp)/P MOD 1.)

; Eccentric anomaly.
EE=keplerec2(M,e)
;print,'EE=',EE

; True anomaly.
theta=2.*atan(sqrt((1.+e)/(1.-e))*tan(EE/2.))
;print,'theta=',theta

; Radial velocity
; rv = gamma + K*[ cos(theta(t)+omega) + e*cos(omega) ]
rv=gamma+K*(cos(theta+omega)+e*cos(omega))

return,rv
end

;-------------------------------------------------------------------------------

; This is the chi^2 function for a single-lined binary.
function chisqSL,param
common datos,t1,rv1,srv1,t2,rv2,srv2

rvcalc1=radvel2(t1,param)

return,total(((rvcalc1-rv1)/srv1)^2)
end




; This is the chi^2 function for a double-lined binary.
function chisqDL,param
common datos,t1,rv1,srv1,t2,rv2,srv2

param1=param[[0,1,2,3,4,5]]
rvcalc1=radvel2(t1,param1)

param2=param[[0,1,2,3,4,6]]
param2[3]=(param2[3]+180d) mod 360.	; omega+180º for the secondary.
rvcalc2=radvel2(t2,param2)

return,total(((rvcalc1-rv1)/srv1)^2) + total(((rvcalc2-rv2)/srv2)^2)

end

;-------------------------------------------------------------------------------

function ConfIntShort,h,theta,C,theta_max,maxhisto=maxhisto
; Function to compute the shortes confidence interval from the histogram.
; h = histogram.
; theta = position of each bin in the histogram.
; C = confidence interval in percentage.
; theta_max = histogram maximum or expected value.
; /maxhisto = locate the bin wich contain the maximum of the histogram.

nbins=n_elements(h)

; Locate the maximum of the histogram.
if keyword_set(maxhisto) then begin
  kk=max(h,ind)
  theta_max=theta(ind)
endif


; CDF from the histogram.
cdf=TOTAL(h, /CUMULATIVE)/TOTAL(FLOAT(h))

; Compute the shortest interval.
eps=0.01*abs(theta[1]-theta[0])	; un trozo do binsize.
dmin=abs(theta[0]-theta[nbins-1]); separación inicial = dominio do histograma.
theta_inf_test=theta[0]

kk=0
while 1 do begin
  ; Add one step.
  theta_inf_test=theta_inf_test+eps
  
  ; Compute the lower cdf of the confidence interval.
  cdf_inf_test=interpol(cdf,theta,theta_inf_test)
  
  ; The diference between the upper cdf and the lower cdf always is C.
  cdf_sup_test=cdf_inf_test+C
  ;print,'theta_inf_test=',theta_inf_test,' cdf_inf_test=',cdf_inf_test,' cdf_sup_test=',cdf_sup_test  

  if cdf_sup_test ge 1.0 then begin
  ; When the upper computed cdf is greater than 1 then end the loop.
    break
  endif else begin
    ; Compute the theta associated to the upper cdf and the width of the confidence interval.
    theta_sup_test=interpol(theta,cdf,cdf_sup_test)
    d=abs(theta_sup_test-theta_inf_test)	; Compute the difference.
    
    ; When the difference is smaller then save the data.
    if d lt dmin then begin
      dmin=d
      theta_inf=theta_inf_test
      theta_sup=theta_sup_test
      cdf_inf=cdf_inf_test
      cdf_sup=cdf_sup_test
    endif
  endelse

endwhile

porcent=abs(cdf_inf-cdf_sup)
return,[theta_inf,theta_max,theta_sup,porcent]
end

