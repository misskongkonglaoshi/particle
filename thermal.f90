function cal_delh(b)
    !  to calculate the enthalpy difference between air and wall    
        use general
        use grid
        use var0
        use uni0
        use uni1
        use fuel
        implicit none
        integer i,j
        real*8 ttav,ma
        real*8 enthalpy,tmp1,tmp2,tstd,enthalpytotal,enthalpytotal2
        real*8 b
        b=a1(ii)*aref/height!���ڴ�����
        do j=1,inj_point!��ȫ���쳤��lmix
            if(phi(j).gt.1.0) then!�����ȴ���1ʱ   cm����ϵ��
                lmix(j)=3.333*cm*dexp(-1.204*phi(j))*b
            else
                lmix(j)=0.179*cm*dexp(1.72*phi(j))*b
            endif
            i=ii
            !��ͬ����ģ��   ˳�� ���� ֧��
            if(inj_mode(j)==1) then
                eta(j)=(x(i)-xinj(j))*lref/lmix(j)
            elseif(inj_mode(j)==2)then
                eta(j)=((x(i)-xinj(j))*lref/lmix(j)+1./(50.+1000.*alf))**alf
            elseif(inj_mode(j)==3)then
                eta(j)=(1.-dexp(-amix*((x(i)-xinj(j))*lref)/lmix(j)))/(1.-dexp(-amix))
            endif
            if(eta(j)>1.d0) eta(j)=1.d0
        enddo
        ttav=0.5d0*(t(1)+t(0))
        ma=(u(1)+u(0))/(c(1)+c(0))
        ttav=ttav*(1.+(shr(1)-1.)/2.0*ma*ma)
        tmp1=0.
        tstd=298.
        do j=1,inj_point
            tmp1=tmp1+phi(j)*stoi(j)*(enthalpy(tinj(j)*tref,fuel_type(j)+1)-enthalpy(tstd,fuel_type(j)+1))
        enddo
        tmp2=0.
        do j=1,inj_point
            tmp2=tmp2+0.5*eu(j)*eta(j)*phi(j)*stoi(j)
        enddo
        
        cal_delh=enthalpy(ttav*tref,1)+tmp1+tmp2-enthalpy(tw*tref,1)
        
        if (inj_lib(1)==1) then
            cal_delh=enthalpytotal(ii-1,ttav*tref)-enthalpytotal(ii-1,tw*tref)
        endif
        return
        end



function enthalpy(t,mode)  !������ֵ
    implicit none
    integer i,j,mode
    real*8 aw(6),ro,qh(6,6),qh1(6,6),hn2,ho2,hc2h4,hh2,hh2o,hc12h23,t2,t3,t4,t5,t,aw2(8),qh9g1(8,8),qh9g2(8,8),tln
    real*8 hco2,hc10h16,hal,hb,hal2o3,hb2o3,hhbo2,hch4,hals,hbs,hal2o3s,hb2o3s,hhbo2s
    real*8 qhals(8),qhall(8),qhbs(8),qhbl(8),qhal2o3s1(8),qhal2o3s2(8),qhal2o3l(8),qhb2o3s(8),qhb2o3l(8),qhhbo2s(8),qhhbo2l(8)
    data aw/28.014d0,31.998d0,28.032d0,167.0,2.01594d0,18.01534d0/
    data aw2/44.0095d0,136.234d0,26.9815d0,10.811d0,101.961d0,69.62d0,43.8177d0,16.042d0/
    data ((qh1(i,j),j=1,6),i=1,6)/                                     &
        0.02926640e+02,0.14879768e-02,-0.05684760e-05,0.10097038e-09,  &
        -0.06753351e-13,-0.09227977e+04,                               &
        0.03697578e+02,0.06135197e-02,-0.12588420e-06,0.01775281e-09,  &
        -0.11364354e-14,-0.12339301e+04,                               &
        0.03528418e+02,0.11485185e-01,-0.04418385e-04,0.07844600e-08,  &
        -0.05266848e-12,0.04428288e+05,                                &
        2.48802010e+01,7.82500480e-02,-3.15509730e-05,5.78789000e-09,  &
        -3.98279680e-13,-4.31106840e+04,                               &
        0.02991423e+02,0.07000644e-02,-0.05633828e-06,-0.09231578e-10, &
        0.15827519e-14,-0.08350340e+04 ,                               &
        0.02672145e+02,0.03056293e-01,-0.08730260e-05,0.12009964e-09,  &
        -0.06391618e-13,-0.02989921e+06/
    data ((qh(i,j),j=1,6),i=1,6)/                                      &
        0.03298677e+02,0.14082404e-02,-0.03963222e-04,0.05641515e-07,  &
        -0.02444854e-10,-0.10208999e+04,                               &
        0.03212936e+02,0.11274864e-02,-0.05756150e-05,0.13138773e-08,  &
        -0.08768554e-11,-0.10052490e+04,                               &
        -0.08614880e+01,0.02796162e+00,-0.03388677e-03,0.02785152e-06, &
        -0.09737879e-10,0.05573046e+05,                                &
        2.08692170e+00,1.33149650e-01,-8.11574520e-05,2.94092860e-08,  &
        -6.51952130e-12,-3.59128140e+04,                               &
        0.03298124e+02,0.08249441e-02,-0.08143015e-05,-0.09475434e-09, &
        0.04134872e-11,-0.10125209e+04,                                &
        0.03386842e+02,0.03474982e-01,-0.06354696e-04,0.06968581e-07,  &
        -0.02506588e-10,-0.03020811e+06/
    data ((qh9g1(i,j),j=1,8),i=1,8)/                                                        &
        4.943650540e+04,-6.264116010e+02, 5.301725240e+00, 2.503813816e-03,-2.127308728e-07,&
        -7.689988780e-10, 2.849677801e-13,-4.528198460e+04,                                 &
        -7.310769440e+05, 1.521764245e+04,-1.139312644e+02,4.281501620e-01,-5.218740440e-04,&
        3.357233400e-07,-8.805750980e-11 ,-8.067482120e+04,                                 &
        -2.920820938e+04,1.167751876e+02,2.356906505e+00,7.737231520e-05,-1.529455262e-08,  &
        -9.971670260e-13, 5.053278264e-16,3.823288650e+04,                                  &
        -1.072659610e+05,3.225307160e+02,2.126407232e+00,2.106579339e-04,-5.937129160e-08,  &
        7.377427990e-12,-2.282443381e-16,6.643413100e+04,                                   &
        -7.443374320e+03,8.829004210e+01, 5.264662640e+00, 2.507678848e-02,-3.434541650e-05,&
        2.302516980e-08,-6.122529280e-12,-6.872685950e+04,                                   &
        7.379611910e+04,-1.263620592e+03, 1.072681512e+01, 3.841383720e-04, 5.976058380e-06,&
        -6.552891350e-09, 2.123951064e-12,-9.628183140e+04,                                 &
        6.225087470e+03,7.566153690e+01,1.253406833e+00,1.748006535e-02,-1.982688351e-05,   &
        1.229656460e-08,-3.153609847e-12,-6.878588780e+04,                                  &
        -1.766850998e+05,2.786181020e+03,-1.202577850e+01,3.917619290e-02,-3.619054430e-05, &
        2.026853043e-08,-4.976705490e-12,-2.331314360e+04/
        
    data ((qh9g2(i,j),j=1,8),i=1,8)/                                                        &
        1.176962419e+05,-1.788791477e+03, 8.291523190e+00,-9.223156780e-05, 4.863676880e-09,&
        -1.891053312e-12, 6.330036590e-16,-3.908350590e+04,                                  &
        1.220329594e+07,-5.794846240e+04 ,1.092281156e+02,-1.082406215e-02,2.034992622e-06, &
        -2.052060369e-10, 8.575760210e-15 , 3.257334050e+05,                                &
        -2.920820938e+04,1.167751876e+02,2.356906505e+00,7.737231520e-05,-1.529455262e-08,  &
        -9.971670260e-13, 5.053278264e-16,3.823288650e+04,                                  &
        -1.072659610e+05,3.225307160e+02,2.126407232e+00,2.106579339e-04,-5.937129160e-08,  &
        7.377427990e-12,-2.282443381e-16,6.643413100e+04,                                   &
        -2.777784969e+05,-4.917465930e+02,1.386703888e+01,-1.469381940e-04,3.250406490e-08, &
        -3.730867350e-12,1.730444284e-16,-6.790757850e+04,                                  &
        3.905035300e+05,-3.691348210e+03,1.555502598e+01,-9.707645510e-04, 2.068887872e-07, &
        -2.310858356e-11,1.050136734e-15,-8.263054410e+04,                                  &
        1.049369185e+06,-4.479145480e+03,1.197755861e+01,-4.735743400e-04, 6.080207140e-08, &
        -3.641565440e-12, 6.155973170e-17,-4.221149470e+04,                                 &
        3.730042760e+06,-1.383501485e+04, 2.049107091e+01,-1.961974759e-03,4.727313040e-07, &
        -3.728814690e-11, 1.623737207e-15,7.532066910e+04/
    qhals=[-6.251811430e+04, 6.343934350e+02,-7.131883820e-01, 1.088725280e-02,-1.458741820e-05,&
        9.961160880e-09,-1.774928010e-12,-3.985439320e+03]!200-934
    qhall=[0.000000000D+00, 0.000000000D+00, 3.818625510D+00, 0.000000000D+00, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-9.576323160D+01]!934-6000
    qhbs=[-8.697700220D+02,-8.050405960D+02, 4.079712880D+00,-6.423381350D-04, 4.846017800D-07,&
        -1.252780673D-10, 1.335923595D-14,3.397919930D+03]!600-2350
    qhbl=[ 0.000000000D+00, 0.000000000D+00, 3.818625511D+00, 0.000000000D+00, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,3.360603140D+03]!2350-6000
    qhal2o3s1=[-6.042087868D+05, 0.000000000D+00, 1.475480816D+01, 8.272285438D-04, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-2.079235447D+05]!500-1200
    qhal2o3s2=[0.000000000D+00,0.000000000D+00, 1.293774378D+01, 1.992781294D-03, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-2.060787581D+05]!1200-2327
    qhal2o3l=[0.000000000D+00, 0.000000000D+00, 1.959225499D+01, 0.000000000D+00, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-2.027701571D+05]!2327-6000
    qhb2o3s=[-5.595297380D+04, 1.311214190D+03,-1.178535942D+01, 7.702795250D-02,-9.740126500D-05,&
        4.692119720D-08, 1.804813810D-12, -1.599672923D+05]!100-723
    qhb2o3l=[3.774124994D+05, 0.000000000D+00, 1.528015481D+01, 0.000000000D+00, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-1.562115789D+05]!723-6000
    qhhbo2s=[0.000000000D+00, 0.000000000D+00, 4.478916978D+00, 7.045514408D-03, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-9.841912440D+04]!300-509
    qhhbo2l=[0.000000000D+00, 0.000000000D+00, 1.262852531D+01, 0.000000000D+00, 0.000000000D+00,&
        0.000000000D+00, 0.000000000D+00,-9.993471070D+04]!509-6000
    ro=8314.d0
    tln=log(t)
    t2=t**2
    t3=t**3
    t4=t**4
    t5=t**5
    if(t.lt.1000.) then
        hn2=ro/aw(1)*(qh(1,1)*t+qh(1,2)*t2/2.+qh(1,3)*t3/3.+qh(1,4)*t4/4.+qh(1,5)*t5/5.+qh(1,6))
        ho2=ro/aw(2)*(qh(2,1)*t+qh(2,2)*t2/2.+qh(2,3)*t3/3.+qh(2,4)*t4/4.+qh(2,5)*t5/5.+qh(2,6))
        hc2h4=ro/aw(3)*(qh(3,1)*t+qh(3,2)*t2/2.+qh(3,3)*t3/3.+qh(3,4)*t4/4.+qh(3,5)*t5/5.+qh(3,6))
        hc12h23=ro/aw(4)*(qh(4,1)*t+qh(4,2)*t2/2.+qh(4,3)*t3/3.+qh(4,4)*t4/4.+qh(4,5)*t5/5.+qh(4,6))
        hh2=ro/aw(5)*(qh(5,1)*t+qh(5,2)*t2/2.+qh(5,3)*t3/3.+qh(5,4)*t4/4.+qh(5,5)*t5/5.+qh(5,6))
        hh2o=ro/aw(6)*(qh(6,1)*t+qh(6,2)*t2/2.+qh(6,3)*t3/3.+qh(6,4)*t4/4.+qh(6,5)*t5/5.+qh(6,6))
    else
        hn2=ro/aw(1)*(qh1(1,1)*t+qh1(1,2)*t2/2.+qh1(1,3)*t3/3.+qh1(1,4)*t4/4.+qh1(1,5)*t5/5.+qh1(1,6))
        ho2=ro/aw(2)*(qh1(2,1)*t+qh1(2,2)*t2/2.+qh1(2,3)*t3/3.+qh1(2,4)*t4/4.+qh1(2,5)*t5/5.+qh1(2,6))
        hc2h4=ro/aw(3)*(qh1(3,1)*t+qh1(3,2)*t2/2.+qh1(3,3)*t3/3.+qh1(3,4)*t4/4.+qh1(3,5)*t5/5.+qh1(3,6))
        hc12h23=ro/aw(4)*(qh1(4,1)*t+qh1(4,2)*t2/2.+qh1(4,3)*t3/3.+qh1(4,4)*t4/4.+qh1(4,5)*t5/5.+qh1(4,6))
        hh2=ro/aw(5)*(qh1(5,1)*t+qh1(5,2)*t2/2.+qh1(5,3)*t3/3.+qh1(5,4)*t4/4.+qh1(5,5)*t5/5.+qh1(5,6))
        hh2o=ro/aw(6)*(qh1(6,1)*t+qh1(6,2)*t2/2.+qh1(6,3)*t3/3.+qh1(6,4)*t4/4.+qh1(6,5)*t5/5.+qh1(6,6))
    endif
    if(t.lt.1000.) then
        hco2=ro/aw2(1)*(-qh9g1(1,1)/t+qh9g1(1,2)*tln+qh9g1(1,3)*t+qh9g1(1,4)*t2/2+qh9g1(1,5)*t3/3+qh9g1(1,6)*t4/4+qh9g1(1,7)*t5/5+qh9g1(1,8))
        hc10h16=ro/aw2(2)*(-qh9g1(2,1)/t+qh9g1(2,2)*tln+qh9g1(2,3)*t+qh9g1(2,4)*t2/2+qh9g1(2,5)*t3/3+qh9g1(2,6)*t4/4+qh9g1(2,7)*t5/5+qh9g1(2,8))
        hal=ro/aw2(3)*(-qh9g1(3,1)/t+qh9g1(3,2)*tln+qh9g1(3,3)*t+qh9g1(3,4)*t2/2+qh9g1(3,5)*t3/3+qh9g1(3,6)*t4/4+qh9g1(3,7)*t5/5+qh9g1(3,8))
        hb=ro/aw2(4)*(-qh9g1(4,1)/t+qh9g1(4,2)*tln+qh9g1(4,3)*t+qh9g1(4,4)*t2/2+qh9g1(4,5)*t3/3+qh9g1(4,6)*t4/4+qh9g1(4,7)*t5/5+qh9g1(4,8))
        hal2o3=ro/aw2(5)*(-qh9g1(5,1)/t+qh9g1(5,2)*tln+qh9g1(5,3)*t+qh9g1(5,4)*t2/2+qh9g1(5,5)*t3/3+qh9g1(5,6)*t4/4+qh9g1(5,7)*t5/5+qh9g1(5,8))
        hb2o3=ro/aw2(6)*(-qh9g1(6,1)/t+qh9g1(6,2)*tln+qh9g1(6,3)*t+qh9g1(6,4)*t2/2+qh9g1(6,5)*t3/3+qh9g1(6,6)*t4/4+qh9g1(6,7)*t5/5+qh9g1(6,8))
        hhbo2=ro/aw2(7)*(-qh9g1(7,1)/t+qh9g1(7,2)*tln+qh9g1(7,3)*t+qh9g1(7,4)*t2/2+qh9g1(7,5)*t3/3+qh9g1(7,6)*t4/4+qh9g1(7,7)*t5/5+qh9g1(7,8))
        hch4=ro/aw2(8)*(-qh9g1(8,1)/t+qh9g1(8,2)*tln+qh9g1(8,3)*t+qh9g1(8,4)*t2/2+qh9g1(8,5)*t3/3+qh9g1(8,6)*t4/4+qh9g1(8,7)*t5/5+qh9g1(8,8))
    else
        hco2=ro/aw2(1)*(-qh9g2(1,1)/t+qh9g2(1,2)*tln+qh9g2(1,3)*t+qh9g2(1,4)*t2/2+qh9g2(1,5)*t3/3+qh9g2(1,6)*t4/4+qh9g2(1,7)*t5/5+qh9g2(1,8))
        hc10h16=ro/aw2(2)*(-qh9g2(2,1)/t+qh9g2(2,2)*tln+qh9g2(2,3)*t+qh9g2(2,4)*t2/2+qh9g2(2,5)*t3/3+qh9g2(2,6)*t4/4+qh9g2(2,7)*t5/5+qh9g2(2,8))
        hal=ro/aw2(3)*(-qh9g2(3,1)/t+qh9g2(3,2)*tln+qh9g2(3,3)*t+qh9g2(3,4)*t2/2+qh9g2(3,5)*t3/3+qh9g2(3,6)*t4/4+qh9g2(3,7)*t5/5+qh9g2(3,8))
        hb=ro/aw2(4)*(-qh9g2(4,1)/t+qh9g2(4,2)*tln+qh9g2(4,3)*t+qh9g2(4,4)*t2/2+qh9g2(4,5)*t3/3+qh9g2(4,6)*t4/4+qh9g2(4,7)*t5/5+qh9g2(4,8))
        hal2o3=ro/aw2(5)*(-qh9g2(5,1)/t+qh9g2(5,2)*tln+qh9g2(5,3)*t+qh9g2(5,4)*t2/2+qh9g2(5,5)*t3/3+qh9g2(5,6)*t4/4+qh9g2(5,7)*t5/5+qh9g2(5,8))
        hb2o3=ro/aw2(6)*(-qh9g2(6,1)/t+qh9g2(6,2)*tln+qh9g2(6,3)*t+qh9g2(6,4)*t2/2+qh9g2(6,5)*t3/3+qh9g2(6,6)*t4/4+qh9g2(6,7)*t5/5+qh9g2(6,8))
        hhbo2=ro/aw2(7)*(-qh9g2(7,1)/t+qh9g2(7,2)*tln+qh9g2(7,3)*t+qh9g2(7,4)*t2/2+qh9g2(7,5)*t3/3+qh9g2(7,6)*t4/4+qh9g2(7,7)*t5/5+qh9g2(7,8))
        hch4=ro/aw2(8)*(-qh9g2(8,1)/t+qh9g2(8,2)*tln+qh9g2(8,3)*t+qh9g2(8,4)*t2/2+qh9g2(8,5)*t3/3+qh9g2(8,6)*t4/4+qh9g2(8,7)*t5/5+qh9g2(8,8))
    endif
    select case(mode)
    case(1)!air
        enthalpy=0.233*ho2+0.767*hn2
    case(2)!yixi
        enthalpy=hc2h4
    case(3)!����
        enthalpy=hh2
    case (4)!ú��
        enthalpy=hc12h23
    case(5)!ˮ
        enthalpy=hh2o
    case(6)!����
        enthalpy=hn2
    case(7)!����
        enthalpy=ho2
    case(8)!������̼
        enthalpy=hco2
    case(9)!jp-10
        enthalpy=hc10h16
    case(10)!Al
        enthalpy=hal
    case(11)!B
        enthalpy=hb
    case(12)!al2o3
        enthalpy=hal2o3
    case(13)!b2o3
        enthalpy=hb2o3
    case(14)!hbo2
        enthalpy=hhbo2
    case(15)!ch4
        enthalpy=hch4
    case(16)!al-s
        if(t.lt.934.) then
            hals=ro/aw2(3)*(-qhals(1)/t+qhals(2)*tln+qhals(3)*t+qhals(4)*t2/2+qhals(5)*t3/3+qhals(6)*t4/4+qhals(7)*t5/5+qhals(8))
        else
            hals=ro/aw2(3)*(-qhall(1)/t+qhall(2)*tln+qhall(3)*t+qhall(4)*t2/2+qhall(5)*t3/3+qhall(6)*t4/4+qhall(7)*t5/5+qhall(8))
        endif        
        enthalpy=hals
    case(17)!b-s
        if(t.lt.2250.) then
            hbs=ro/aw2(4)*(-qhbs(1)/t+qhbs(2)*tln+qhbs(3)*t+qhbs(4)*t2/2+qhbs(5)*t3/3+qhbs(6)*t4/4+qhbs(7)*t5/5+qhbs(8))
        elseif(t.lt.2350.) then
            hbs=(2350-t)/100*ro/aw2(4)*(-qhbs(1)/t+qhbs(2)*tln+qhbs(3)*t+qhbs(4)*t2/2+qhbs(5)*t3/3+qhbs(6)*t4/4+qhbs(7)*t5/5+qhbs(8))+(t-2250)/100*ro/aw2(4)*(-qhbl(1)/t+qhbl(2)*tln+qhbl(3)*t+qhbl(4)*t2/2+qhbl(5)*t3/3+qhbl(6)*t4/4+qhbl(7)*t5/5+qhbl(8))
        else
            hbs=ro/aw2(4)*(-qhbl(1)/t+qhbl(2)*tln+qhbl(3)*t+qhbl(4)*t2/2+qhbl(5)*t3/3+qhbl(6)*t4/4+qhbl(7)*t5/5+qhbl(8))
        endif        
        enthalpy=hbs 
    case(18)!al2o3-s
        if(t.lt.2427.) then
            if(t.lt.1200.) then
                hal2o3s=ro/aw2(5)*(-qhal2o3s1(1)/t+qhal2o3s1(2)*tln+qhal2o3s1(3)*t+qhal2o3s1(4)*t2/2+qhal2o3s1(5)*t3/3+qhal2o3s1(6)*t4/4+qhal2o3s1(7)*t5/5+qhal2o3s1(8))
            elseif(t.lt.2227.) then
                hal2o3s=ro/aw2(5)*(-qhal2o3s2(1)/t+qhal2o3s2(2)*tln+qhal2o3s2(3)*t+qhal2o3s2(4)*t2/2+qhal2o3s2(5)*t3/3+qhal2o3s2(6)*t4/4+qhal2o3s2(7)*t5/5+qhal2o3s2(8))
            else
                hal2o3s=(2427.d0-t)/200.d0*(ro/aw2(5)*(-qhal2o3s2(1)/t+qhal2o3s2(2)*tln+qhal2o3s2(3)*t+qhal2o3s2(4)*t2/2+qhal2o3s2(5)*t3/3+qhal2o3s2(6)*t4/4+qhal2o3s2(7)*t5/5+qhal2o3s2(8)))+(t-2227.d0)/200.d0*ro/aw2(5)*(-qhal2o3l(1)/t+qhal2o3l(2)*tln+qhal2o3l(3)*t+qhal2o3l(4)*t2/2+qhal2o3l(5)*t3/3+qhal2o3l(6)*t4/4+qhal2o3l(7)*t5/5+qhal2o3l(8))
            endif
        else
            hal2o3s=ro/aw2(5)*(-qhal2o3l(1)/t+qhal2o3l(2)*tln+qhal2o3l(3)*t+qhal2o3l(4)*t2/2+qhal2o3l(5)*t3/3+qhal2o3l(6)*t4/4+qhal2o3l(7)*t5/5+qhal2o3l(8))
        endif        
        enthalpy=hal2o3s
    case(19)!b2o3-s
        if(t.lt.723.) then
            hb2o3s=ro/aw2(6)*(-qhb2o3s(1)/t+qhb2o3s(2)*tln+qhb2o3s(3)*t+qhb2o3s(4)*t2/2+qhb2o3s(5)*t3/3+qhb2o3s(6)*t4/4+qhb2o3s(7)*t5/5+qhb2o3s(8))  
        else
            hb2o3s=ro/aw2(6)*(-qhb2o3l(1)/t+qhb2o3l(2)*tln+qhb2o3l(3)*t+qhb2o3l(4)*t2/2+qhb2o3l(5)*t3/3+qhb2o3l(6)*t4/4+qhb2o3l(7)*t5/5+qhb2o3l(8))
        endif        
        enthalpy=hb2o3s
    case(20)!HBO2-s
        if(t.lt.509.) then
            hhbo2s=ro/aw2(7)*(-qhhbo2s(1)/t+qhhbo2s(2)*tln+qhhbo2s(3)*t+qhhbo2s(4)*t2/2+qhhbo2s(5)*t3/3+qhhbo2s(6)*t4/4+qhhbo2s(7)*t5/5+qhhbo2s(8))  
        else
            hhbo2s=ro/aw2(7)*(-qhhbo2l(1)/t+qhhbo2l(2)*tln+qhhbo2l(3)*t+qhhbo2l(4)*t2/2+qhhbo2l(5)*t3/3+qhhbo2l(6)*t4/4+qhhbo2l(7)*t5/5+qhhbo2l(8))
        endif        
        enthalpy=hhbo2s
    case(21)!
        case default
        !print*,"no such material'enthalpy"
    end select
    end function enthalpy
    
    !ȷ����ֵ���ֵ����   ��������*��


function enthalpytotal(i,tt)
    use var0
    implicit none
    integer i
    real*8 tt,enthalpy,enthalpytotal
    enthalpytotal=mgo2(i)*enthalpy(tt,7) + mgn2(i)*enthalpy(tt,6) + mgco2(i)*enthalpy(tt,8) + mgh2o(i)*enthalpy(tt,5) + mgjp10(i)*enthalpy(tt,9) + mgal2o3(i)*enthalpy(tt,12) + mgb2o3(i)*enthalpy(tt,13) + mghbo2(i)*enthalpy(tt,14)
    end function enthalpytotal
    
function enthalpytotal2(i,tt)
    use var0
    implicit none
    integer i
    real*8 tt,enthalpy,enthalpytotal2
    enthalpytotal2=msal(i)*enthalpy(tt,16) + msb(i)*enthalpy(tt,17) + msal2o3(i)*enthalpy(tt,18) + msb2o3(i)*enthalpy(tt,19) + mshbo2(i)*enthalpy(tt,20)
    end function enthalpytotal2
    
function cp_cal(t,mode)!���㲻ͬ��֣�mode��cp   
    implicit none
    integer mode,i,j
    real*8::qh1(7,6),qh(7,6)
    real*8::cp_cal,t,aw(7),rgl,r
    data aw/28.014d0,31.998d0,28.032d0,167.0,2.01594d0,18.01534d0,28.032d0/      !n2,o2,c2h4,c12h23,h2,h2o,co2
    data ((qh1(i,j),j=1,6),i=1,7)/&
        0.02926640e+02,0.14879768e-02,-0.05684760e-05,0.10097038e-09,-0.06753351e-13,-0.09227977e+04,&
        0.03697578e+02,0.06135197e-02,-0.12588420e-06,0.01775281e-09,-0.11364354e-14,-0.12339301e+04,&
        0.03528418e+02,0.11485185e-01,-0.04418385e-04,0.07844600e-08,-0.05266848e-12,0.044282880e+05,&
        2.48802010e+01,7.82500480e-02,-3.15509730e-05,5.78789000e-09,-3.98279680e-13,-4.31106840e+04,&
        0.02991423e+02,0.07000644e-02,-0.05633828e-06,-0.09231578e-10,0.15827519e-14,-0.08350340e+04,&
        0.02672145e+02,0.03056293e-01,-0.08730260e-05,0.12009964e-09,-0.06391618e-13,-0.02989921e+06,&
        0.04453623E+02,0.03140168E-01,-0.12784105E-05,0.02393996E-08,-0.16690333E-13,-0.04896696E+06/
    data ((qh(i,j),j=1,6),i=1,7)/&
        0.03298677e+02,0.14082404e-02,-0.03963222e-04,0.05641515e-07,-0.02444854e-10,-0.10208999e+04,&
        0.03212936e+02,0.11274864e-02,-0.05756150e-05,0.13138773e-08,-0.08768554e-11,-0.10052490e+04,&
        -0.08614880e+01,0.02796162e+00,-0.03388677e-03,0.02785152e-06,-0.09737879e-10,0.05573046e+05,&
        2.08692170e+00,1.33149650e-01,-8.11574520e-05,2.94092860e-08,-6.51952130e-12,-3.59128140e+04,&
        0.03298124e+02,0.08249441e-02,-0.08143015e-05,-0.09475434e-09,0.04134872e-11,-0.10125209e+04,&
        0.03386842e+02,0.03474982e-01,-0.06354696e-04,0.06968581e-07,-0.02506588e-10,-0.03020811e+06,&
        0.04453623E+02,0.03140168E-01,-0.12784105E-05,0.02393996E-08,-0.16690333E-13,-0.04896696E+06/
    r=8314.d0
    rgl=8314.d0/aw(mode)
    if (t<=1000.d0) then
        cp_cal=rgl*(qh(mode,1)+qh(mode,2)*t+qh(mode,3)*t*t+qh(mode,4)*t**3+qh(mode,5)*t**4)
    else
        cp_cal=rgl*(qh1(mode,1)+qh1(mode,2)*t+qh1(mode,3)*t*t+qh1(mode,4)*t**3+qh1(mode,5)*t**4)
    endif
    end function cp_cal
    
function cp_cal2(t,mode)!��ú�� �� �� ��cp����
    implicit none
    integer mode,i,j
    real*8 aw(8),ro,qh(8,8),qh2(8,8),t2,t3,t4,t5,t,tln,cp_cal2
    data aw/31.998d0,28.014d0,44.0095d0,18.01534d0,136.234d0,69.62d0,43.8177d0,101.961d0/!o2 n2 co2 h2o jp10 b2o3 hbo2 al2o3
    data ((qh(i,j),j=1,8),i=1,8)/                                                           &
        -3.425563420e+04,4.847000970e+02,1.119010961e+00,4.293889240e-03,-6.836300520e-07,  &
        -2.023372700e-09, 1.039040018e-12,-3.391454870e+03,                                 &
        2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,  &
        -9.625793620e-09, 2.519705809e-12,7.108460860e+02,                                  &
        4.943650540e+04,-6.264116010e+02, 5.301725240e+00, 2.503813816e-03,-2.127308728e-07,&
        -7.689988780e-10, 2.849677801e-13,-4.528198460e+04,                                 &
        -3.947960830D+04,5.755731020D+02,9.317826530D-01,7.222712860D-03,-7.342557370D-06,  &
        4.955043490D-09,-1.336933246D-12,-3.303974310D+04,                                  &
        -7.310769440e+05, 1.521764245e+04,-1.139312644e+02,4.281501620e-01,-5.218740440e-04,&
        3.357233400e-07,-8.805750980e-11 ,-8.067482120e+04,                                 &
        7.379611910e+04,-1.263620592e+03, 1.072681512e+01, 3.841383720e-04, 5.976058380e-06,&
        -6.552891350e-09, 2.123951064e-12,-9.628183140e+04,                                 &
        6.225087470e+03,7.566153690e+01,1.253406833e+00,1.748006535e-02,-1.982688351e-05,   &
        1.229656460e-08,-3.153609847e-12,-6.878588780e+04,                                  &
        -7.443374320e+03,8.829004210e+01,5.264662640e+00,2.507678848e-02,-3.434541650e-05,  &
        2.302516980e-08,-6.122529280e-12,-6.872685950D+04/
    data ((qh2(i,j),j=1,8),i=1,8)/                                                        &
        -1.037939022e+06,2.344830282e+03,1.819732036e+00,1.267847582e-03,-2.188067988e-07,  &
        2.053719572e-11,-8.193467050e-16,-1.689010929e+04,                                  &
        5.877124060D+05,-2.239249073D+03,6.066949220D+00,-6.139685500D-04,1.491806679D-07,  &
        -1.923105485D-11, 1.061954386D-15,1.283210415D+04,                                  &
        1.176962419e+05,-1.788791477e+03, 8.291523190e+00,-9.223156780e-05, 4.863676880e-09,&
        -1.891053312e-12, 6.330036590e-16,-3.908350590e+04,                                 &
        1.034972096D+06,-2.412698562D+03,4.646110780D+00,2.291998307D-03,-6.836830480D-07,  &
        9.426468930D-11,-4.822380530D-15,-1.384286509D+04,                                  &
        1.220329594e+07,-5.794846240e+04 ,1.092281156e+02,-1.082406215e-02,2.034992622e-06, &
        -2.052060369e-10, 8.575760210e-15 , 3.257334050e+05,                                &
        3.905035300e+05,-3.691348210e+03,1.555502598e+01,-9.707645510e-04, 2.068887872e-07, &
        -2.310858356e-11,1.050136734e-15,-8.263054410e+04,                                  &
        1.049369185e+06,-4.479145480e+03,1.197755861e+01,-4.735743400e-04, 6.080207140e-08, &
        -3.641565440e-12, 6.155973170e-17,-4.221149470e+04,                                 &
        -2.777784969e+05,-4.917465930e+02,1.386703888e+01,-1.469381940e-04,3.250406490e-08, &
        -3.730867350e-12, 1.730444284e-16,-6.790757850D+04/
    ro=8314.d0
    t2=t**2
    t3=t**3
    t4=t**4
    if (t<=1000.d0) then
        cp_cal2=ro/aw(mode)*(qh(mode,1)/t2+qh(mode,2)/t+qh(mode,3)+qh(mode,4)*t+qh(mode,5)*t2+qh(mode,6)*t3+qh(mode,7)*t4)
    else
        cp_cal2=ro/aw(mode)*(qh2(mode,1)/t2+qh2(mode,2)/t+qh2(mode,3)+qh2(mode,4)*t+qh2(mode,5)*t2+qh2(mode,6)*t3+qh2(mode,7)*t4)
    endif
    end function cp_cal2
  

subroutine heat_ratio!���¼�����ȱ�
    use cal
    use var0
    use uni0
    use uni1
    use combine
    use fuel
    implicit none
    integer i,j
    real*8 cp_cal,cpm,cp_cal2
    real*8::y_o2,y_n2
    y_o2=0.233d0
    y_n2=0.767d0
    do i=-1,ii+1
        cpm=0.d0
        if (inj_lib(1)==1) then
            if (i<1) then!��ע�����������㣿   rflow ƽ��Ħ������
                cpm=cp_cal2(t(i)*tref,1)*mgo2(1) + cp_cal2(t(i)*tref,2)*mgn2(1) + cp_cal2(t(i)*tref,3)*mgco2(1) + cp_cal2(t(i)*tref,4)*mgh2o(1) + cp_cal2(t(i)*tref,5)*mgjp10(1) + cp_cal2(t(i)*tref,8)*mgal2o3(1) + cp_cal2(t(i)*tref,6)*mgb2o3(1) + cp_cal2(t(i)*tref,7)*mghbo2(1)
                rflow(i)=(mgo2(1) + mgn2(1) + mgco2(1) + mgh2o(1) + mgjp10(1) + mgal2o3(1) + mgb2o3(1) + mghbo2(1))/(mgo2(1)/31.998d0+mgn2(1)/28.014d0+mgco2(1)/44.0095d0+mgh2o(1)/18.01534d0+mgjp10(1)/136.234d0+mgal2o3(1)/101.96d0+mgb2o3(1)/69.62d0+mghbo2(1)/43.8177d0)
            elseif(i>ii-1) then!
                rflow(i)=(mgo2(ii-1) + mgn2(ii-1) + mgco2(ii-1) + mgh2o(ii-1) + mgjp10(ii-1) + mgal2o3(ii-1) + mgb2o3(ii-1) + mghbo2(ii-1))/(mgo2(ii-1)/31.998d0+mgn2(ii-1)/28.014d0+mgco2(ii-1)/44.0095d0+mgh2o(ii-1)/18.01534d0+mgjp10(ii-1)/136.234d0+mgal2o3(ii-1)/101.96d0+mgb2o3(ii-1)/69.62d0+mghbo2(ii-1)/43.8177d0)
                cpm=cp_cal2(t(i)*tref,1)*mgo2(ii-1) + cp_cal2(t(i)*tref,2)*mgn2(ii-1) + cp_cal2(t(i)*tref,3)*mgco2(ii-1) + cp_cal2(t(i)*tref,4)*mgh2o(ii-1) + cp_cal2(t(i)*tref,5)*mgjp10(ii-1) + cp_cal2(t(i)*tref,8)*mgal2o3(ii-1) + cp_cal2(t(i)*tref,6)*mgb2o3(ii-1) + cp_cal2(t(i)*tref,7)*mghbo2(ii-1)
            else
                cpm=cp_cal2(t(i)*tref,1)*mgo2(i) + cp_cal2(t(i)*tref,2)*mgn2(i) + cp_cal2(t(i)*tref,3)*mgco2(i) + cp_cal2(t(i)*tref,4)*mgh2o(i) + cp_cal2(t(i)*tref,5)*mgjp10(i) + cp_cal2(t(i)*tref,8)*mgal2o3(i)  + cp_cal2(t(i)*tref,6)*mgb2o3(i) + cp_cal2(t(i)*tref,7)*mghbo2(i)
                rflow(i)=(mgo2(i) + mgn2(i) + mgco2(i) + mgh2o(i) + mgjp10(i) + mgal2o3(i)  + mgb2o3(i) + mghbo2(i))/(mgo2(i)/31.998d0+mgn2(i)/28.014d0+mgco2(i)/44.0095d0+mgh2o(i)/18.01534d0+mgjp10(i)/136.234d0+mgal2o3(i)/101.96d0+mgb2o3(i)/69.62d0+mghbo2(i)/43.8177d0)
            endif
            cp(i)=cpm!ÿ���ڵ��cp
            shr(i)=cp(i)/(cp(i)-8314d0/rflow(i))  !gamma_local
        else
            do j=1,7
                cpm=cpm+mfc(j)*cp_cal(t(i)*tref,j)
            enddo
            cp(i)=cpm
            shr(i)=cp(i)/(cp(i)-r0)
        endif
    enddo
    end subroutine heat_ratio
    
function  heat_cp(i,temp)!���¼������
    use cal
    use var0
    use uni0
    use uni1
    use combine
    use fuel
    implicit none
    integer i,j
    real*8 heat_cp
    real*8 cp_cal,cpm,cp_cal2,temp
    cpm=0.d0
    if (inj_lib(1)==1) then
        cpm=cp_cal2(temp,1)*mgo2(i) + cp_cal2(temp,2)*mgn2(i) + cp_cal2(temp,3)*mgco2(i) + cp_cal2(temp,4)*mgh2o(i) + cp_cal2(temp,5)*mgjp10(i) + cp_cal2(temp,8)*mgal2o3(i) + cp_cal2(temp,6)*mgb2o3(i) + cp_cal2(temp,7)*mghbo2(i)
        heat_cp=cpm
    else
        do j=1,7
            cpm=cpm+mfc(j)*cp_cal(temp,j)
        enddo
        heat_cp=cpm
    endif
    end function heat_cp