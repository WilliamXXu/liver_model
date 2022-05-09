function out=den_exp3(t,y,para2,atp_dec)    %para2=[grwoth_atp_hills_km growth_atp_hills_n  necrosis_k necrosis_atp_hills_km necrosis_atp_hills_n ]  
    densi=y(1)+y(2)+y(3);
    apop=5e-4/3;
    alt_clear=2.43e-4;
    logistic=[0.000499253014622   1.555726168588230 ];

    atp_dec=atp_dec*(atp_dec>0)*(atp_dec<=1)+(atp_dec>1);

    out=[0  %growth, only when (atp_dec<0.5)
    0     %apoptosis
    -para2(1)*HILL(atp_dec,para2(2),para2(3))*densi-para2(5)*HILL(atp_dec,para2(6),para2(7))*densi      %necrosis   
    0
    ];   
     function res=HILL(c,k,n)
        res=1/(1+(k/c)^n);
 end

end