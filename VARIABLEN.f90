module variablen


!! Variablen zur Image-Kontrolle
INTEGER::M,N,O,m_kop,n_kop,o_kop,myrank,numprocs,m_pos,n_pos,o_pos,R,S,T_p&
,m_local,n_local,o_local,image_m_1,image_n_1,image_o_1,image_m_max,image_n_max,&
image_o_max
INTEGER::sync_m_up,sync_m_down,sync_n_up,sync_n_down,sync_o_up,sync_o_down,nlocal
INTEGER,allocatable::sync_m_array(:),sync_n_array(:),sync_o_array(:)

!! Ausgabe-parameter
INTEGER,PARAMETER::prtrows=10
INTEGER,PARAMETER::prtcols=10
LOGICAL::start

!! Felder der Physikgrößen
!! Co-arrays
REAL*8,allocatable::rho(:,:,:)[:,:,:],rho_u(:,:,:)[:,:,:],rho_E(:,:,:)[:,:,:],rho_v(:,:,:)[:,:,:]&
,E(:,:,:)[:,:,:],p(:,:,:)[:,:,:],Epu(:,:,:)[:,:,:],Epv(:,:,:)[:,:,:],Epw(:,:,:)[:,:,:],&
E_k(:,:,:)[:,:,:],vis(:,:,:)[:,:,:],mu(:,:,:)[:,:,:]&
,a(:,:,:)[:,:,:],rho_w(:,:,:)[:,:,:],k(:,:,:)[:,:,:],&
rho_neu(:,:,:)[:,:,:],rho_u_neu(:,:,:)[:,:,:],rho_E_neu(:,:,:)[:,:,:],rho_v_neu(:,:,:)[:,:,:],&
rho_w_neu(:,:,:)[:,:,:]
REAL*8,allocatable::T(:,:,:)[:,:,:],u(:,:,:)[:,:,:],v(:,:,:)[:,:,:],DDX_DDY(:,:,:)[:,:,:],w(:,:,:)[:,:,:]
REAL*8,allocatable::dphi(:,:,:)[:,:,:],solve(:,:,:)[:,:,:],&
dphi1(:,:,:)[:,:,:],dphi2(:,:,:)[:,:,:],dphi3(:,:,:)[:,:,:],dphi4(:,:,:)[:,:,:]
REAL*8,allocatable::du(:,:,:)[:,:,:],dv(:,:,:)[:,:,:]
REAL*8,allocatable::poisson(:,:,:)[:,:,:]
REAL*8,allocatable::uinit(:,:,:)[:,:,:],rinit(:,:,:)[:,:,:],vinit(:,:,:)[:,:,:],tinit(:,:,:)[:,:,:]
REAL*8,allocatable::drudx(:,:,:)[:,:,:],drvdy(:,:,:)[:,:,:],&
dudy(:,:,:)[:,:,:],dvdy(:,:,:)[:,:,:],dudx(:,:,:)[:,:,:],dvdx(:,:,:)[:,:,:],&
rho_bc(:,:,:)[:,:,:],T_bc(:,:,:)[:,:,:],u_bc(:,:,:)[:,:,:],v_bc(:,:,:)[:,:,:]

!! Felder für Gleichungslösung
!in x
REAL*8,allocatable::b_x(:),c_x(:),a_x(:)[:,:,:],f_x(:)[:,:,:],g_x(:)[:,:,:],x_x(:,:,:)[:,:,:],y_x(:,:,:)[:,:,:]
REAL*8,allocatable::a_x_s(:)[:,:,:],f_x_s(:)[:,:,:],g_x_s(:)[:,:,:],y_x_s(:,:,:)[:,:,:],x_x_s(:,:,:)[:,:,:]
REAL*8,allocatable::alpha_x(:,:)[:,:,:],beta_x(:)[:,:,:],z_x(:)
!in y
REAL*8,allocatable::b_y(:),c_y(:),a_y(:)[:,:,:],f_y(:)[:,:,:],g_y(:)[:,:,:],x_y(:,:,:)[:,:,:],y_y(:,:,:)[:,:,:]
REAL*8,allocatable::a_y_s(:)[:,:,:],f_y_s(:)[:,:,:],g_y_s(:)[:,:,:],y_y_s(:,:,:)[:,:,:],x_y_s(:,:,:)[:,:,:]
REAL*8,allocatable::alpha_y(:,:)[:,:,:],beta_y(:)[:,:,:],z_y(:)
!in z
REAL*8,allocatable::b_z(:),c_z(:),a_z(:)[:,:,:],f_z(:)[:,:,:],g_z(:)[:,:,:],x_z(:,:,:)[:,:,:],y_z(:,:,:)[:,:,:]
REAL*8,allocatable::a_z_s(:)[:,:,:],f_z_s(:)[:,:,:],g_z_s(:)[:,:,:],y_z_s(:,:,:)[:,:,:],x_z_s(:,:,:)[:,:,:]
REAL*8,allocatable::alpha_z(:,:)[:,:,:],beta_z(:)[:,:,:],z_z(:)

!! Koeffizienten für kompakte Differenzen
REAL*8::alpha_i,beta_i,gamma_i,alpha_nr,beta_nr,gamma_nr,beta_r,gamma_r
REAL*8::q_i,a_i,b_i,c_i,d_i,e_i,q_nr,b_nr,c_nr,d_nr,q_r,c_r,d_r,e_r
REAL*8::alpha2_i,beta2_i,gamma2_i,alpha2_nr,beta2_nr,gamma2_nr,beta2_r,gamma2_r
REAL*8::q2_i,a2_i,b2_i,c2_i,d2_i,e2_i,q2_nr,b2_nr,c2_nr,d2_nr,q2_r,c2_r,d2_r,e2_r
REAL*8::qxi_i,qxi_nr,qxi_r,qyi_i,qyi_nr,qyi_r,qzi_i,qzi_nr,qzi_r
REAL*8::qxi2_i,qxi2_nr,qxi2_r,qyi2_i,qyi2_nr,qyi2_r,qzi2_i,qzi2_nr,qzi2_r

!! Physik
REAL*8::dx,dy,dz,dt,Ma
REAL*8::R_gas,cv_gas,kappa,cv,cp, T_inf,Pr,gamma
REAL*8::u_0,rho_0,p_0,E_0,vis0
REAL*8::mu0,T0,Suth !mu0, T0 , Sutherland constant
!!
REAL::total_energy[1:*],total_konti[1:*],E_k_total[1:*]
REAL,allocatable::kontinuitaet(:,:)[:]

!!Hilfsvariablen
REAL*8,allocatable::E_plus_p(:,:,:)[:,:,:]

!!spannungstensor
REAL*8,allocatable:: tau_xx(:,:,:)[:,:,:],tau_yy(:,:,:)[:,:,:],tau_zz(:,:,:)[:,:,:]
REAL*8,allocatable:: tau_xy(:,:,:)[:,:,:],tau_xz(:,:,:)[:,:,:],tau_zy(:,:,:)[:,:,:]
REAL*8,allocatable:: delta_u(:,:,:)[:,:,:],delta_v(:,:,:)[:,:,:],delta_w(:,:,:)[:,:,:]
end module
