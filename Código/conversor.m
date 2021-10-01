function res = conversor(val,conversion)
    if strcmp(conversion, 'ua2km')
       res = val*(1.495979e8);
    end
    if strcmp(conversion, 'km2ua')
       res = val/(1.495979e8);
    end
    if strcmp(conversion, 'ua2m')
       res = val*(1.495979e11);
    end
    if strcmp(conversion, 'm2ua')
       res = val/(1.495979e11);
    end
    if strcmp(conversion, 'dia2hora')
        res = val*24;
    end
    if strcmp(conversion, 'hora2min')
        res = val*60;
    end
    if strcmp(conversion, 'min2seg')
        res = val*60;
    end
    if strcmp(conversion, 'dia2seg')
        res = val*86400;
    end  
end