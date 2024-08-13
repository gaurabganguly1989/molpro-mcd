***,cop
memory,16000,M

file,2,cop.wfu,old

set,DKHO=2     ! second-order Douglas-Kroll-Hess Hamiltonian
set,DKHP=1     ! Optimum parametrization (OPT, default)

gprint,basis,orbital,civector;
gthresh, energy  = 1d-9
gthresh, gradient= 1d-6
gthresh, grid    = 1d-9
gthresh, orbital = 1d-8
gthresh, civec   = 1d-8
gexpec,  DM,LOP

symmetry, X Y Z
geomtyp = xyz
angstrom
geometry = cop.doub.xyz

basis = {Co= cc-pVTZ-DK;
         C = cc-pVTZ-DK;
         N = cc-pVTZ-DK;
         H = cc-pVTZ-DK;}

TEXT,"===== STARTING RASSCF CALCULATION ====="
! ^2SA // RAS(15e,16o)
{multi;config,csf;start,2143.2;
      !accuracy, gradient=1d-7, step=1d-4, energy=1d-10;
       !        Ag B3u B2u B1g B1u B2g B3g Au
       frozen,   0,  0,  0,  0,  0,  0,  0, 0;
       closed,  23, 19, 19, 14,  4,  3,  3, 1;
       occ,     27, 19, 19, 16,  7,  6,  6, 2;
      !rotate, 4.5, 5.5, 90;
       wf,      nelec = 187, sym = 6, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state, 1;
       wf,      nelec = 187, sym = 7, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state, 1;
       wf,      nelec = 187, sym = 5, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state,60;
       wf,      nelec = 187, sym = 8, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state,60;
       natorb,2143.2,              cirec = 5143.2;}
       put,molden,nip.cas_16e_16o.3sa.molden;  orbital,2143.2;

TEXT,"===== STARTING MRCI CALCULATION ====="

{mrci; noexc; orbital,2143.2;
       maxiter, 200,200;
       !        Ag B3u B2u B1g B1u B2g B3g Au
       core,    23, 19, 19, 14,  4,  3,  3, 1;
       closed,  23, 19, 19, 14,  4,  3,  3, 1;
       occ,     27, 19, 19, 16,  7,  6,  6, 2;
       wf,      nelec = 187, sym = 6, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state, 1;
       save,6041.2;}

{mrci; noexc; orbital,2143.2;
       maxiter, 200,200;
       !        Ag B3u B2u B1g B1u B2g B3g Au
       core,    23, 19, 19, 14,  4,  3,  3, 1;
       closed,  23, 19, 19, 14,  4,  3,  3, 1;
       occ,     27, 19, 19, 16,  7,  6,  6, 2;
       wf,      nelec = 187, sym = 7, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state, 1;
       save,6051.2;}

{mrci; noexc; orbital,2143.2;
       maxiter, 200,200;
       !        Ag B3u B2u B1g B1u B2g B3g Au
       core,    23, 19, 19, 14,  4,  3,  3, 1;
       closed,  23, 19, 19, 14,  4,  3,  3, 1;
       occ,     27, 19, 19, 16,  7,  6,  6, 2;
       wf,      nelec = 187, sym = 5, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state,60;
       save,6061.2;}

{mrci; noexc; orbital,2143.2;
       maxiter, 200,200;
       !        Ag B3u B2u B1g B1u B2g B3g Au
       core,    23, 19, 19, 14,  4,  3,  3, 1;
       closed,  23, 19, 19, 14,  4,  3,  3, 1;
       occ,     27, 19, 19, 16,  7,  6,  6, 2;
       wf,      nelec = 187, sym = 8, spin = 3, charge = 0; 
       restrict, 0, 2, 26.1, 27.1, 16.4, 6.6, 6.7; 
       restrict,-1,-1, 26.1, 27.1, 16.4, 6.6, 6.7;  
       state,60;
       save,6071.2;}

TEXT,"===== STARTING SO COUPLING CALCULATION ====="

HLSDIAG(001) = -2378.982989464616;   
HLSDIAG(002) = -2378.982989464616;
HLSDIAG(003) = -2378.900456939324;
HLSDIAG(004) = -2378.899305534176;
HLSDIAG(005) = -2378.912900921789;
HLSDIAG(006) = -2378.910106748606;
HLSDIAG(007) = -2378.896352145901;
HLSDIAG(008) = -2378.898400139361;
HLSDIAG(009) = -2378.865695108415;
HLSDIAG(010) = -2378.862927167058;
HLSDIAG(011) = -2378.867698703435;
HLSDIAG(012) = -2378.867495085159;
HLSDIAG(013) = -2378.874538890364;
HLSDIAG(014) = -2378.873133021412;
HLSDIAG(015) = -2378.879645777764;
HLSDIAG(016) = -2378.879353166463;
HLSDIAG(017) = -2378.861216383531;
HLSDIAG(018) = -2378.861062121005;
HLSDIAG(019) = -2378.811967868732;
HLSDIAG(020) = -2378.811786969217;
HLSDIAG(021) = -2378.807720842577;
HLSDIAG(022) = -2378.808130447193;
HLSDIAG(023) = -2378.865806524817;
HLSDIAG(024) = -2378.865025160243;
HLSDIAG(025) = -2378.799782840131;
HLSDIAG(026) = -2378.797760761884;
HLSDIAG(027) = -2378.823660252119;
HLSDIAG(028) = -2378.823289901200;
HLSDIAG(029) = -2378.856139601023;
HLSDIAG(030) = -2378.823029141386;
HLSDIAG(031) = -2378.785907877420;
HLSDIAG(032) = -2378.849997002028;
HLSDIAG(033) = -2378.783312676758;
HLSDIAG(034) = -2378.817309137231;
HLSDIAG(035) = -2378.807096272288;
HLSDIAG(036) = -2378.806523179549;
HLSDIAG(037) = -2378.842271762584;
HLSDIAG(038) = -2378.834472875580;
HLSDIAG(039) = -2378.807995324658;
HLSDIAG(040) = -2378.802316115981;
HLSDIAG(041) = -2378.823455748359;
HLSDIAG(042) = -2378.820884002000;
HLSDIAG(043) = -2378.770816217514;
HLSDIAG(044) = -2378.795055890323;
HLSDIAG(045) = -2378.836373968211;
HLSDIAG(046) = -2378.790188195273;
HLSDIAG(047) = -2378.763606400931;
HLSDIAG(048) = -2378.833830189066;
HLSDIAG(049) = -2378.817242923676;
HLSDIAG(050) = -2378.815880696765;
HLSDIAG(051) = -2378.739897667547;
HLSDIAG(052) = -2378.739648030332;
HLSDIAG(053) = -2378.822209013723;
HLSDIAG(054) = -2378.821894702936;
HLSDIAG(055) = -2378.827972156509;
HLSDIAG(056) = -2378.826591462982;
HLSDIAG(057) = -2378.749212691709;
HLSDIAG(058) = -2378.758332204796;
HLSDIAG(059) = -2378.774432182096;
HLSDIAG(060) = -2378.730217574956;
HLSDIAG(061) = -2378.805169656475;
HLSDIAG(062) = -2378.804099410632;
HLSDIAG(063) = -2378.902061281853;
HLSDIAG(064) = -2378.899145017973;
HLSDIAG(065) = -2378.911549482410;
HLSDIAG(066) = -2378.910224054536;
HLSDIAG(067) = -2378.896843248674;
HLSDIAG(068) = -2378.897909621634;
HLSDIAG(069) = -2378.864658187376;
HLSDIAG(070) = -2378.862162604562;
HLSDIAG(071) = -2378.867958557133;
HLSDIAG(072) = -2378.867564706528;
HLSDIAG(073) = -2378.875258350443;
HLSDIAG(074) = -2378.873899986944;
HLSDIAG(075) = -2378.879544295936;
HLSDIAG(076) = -2378.879449938718;
HLSDIAG(077) = -2378.861122358203;
HLSDIAG(078) = -2378.861151756548;
HLSDIAG(079) = -2378.812278988778;
HLSDIAG(080) = -2378.811746998532;
HLSDIAG(081) = -2378.811248492013;
HLSDIAG(082) = -2378.806544703280;
HLSDIAG(083) = -2378.867232419841;
HLSDIAG(084) = -2378.865424082837;
HLSDIAG(085) = -2378.796354472069;
HLSDIAG(086) = -2378.795330716670;
HLSDIAG(087) = -2378.823416859280;
HLSDIAG(088) = -2378.823319471157;
HLSDIAG(089) = -2378.853889933172;
HLSDIAG(090) = -2378.818034023713;
HLSDIAG(091) = -2378.789027849909;
HLSDIAG(092) = -2378.844327942583;
HLSDIAG(093) = -2378.791255764562;
HLSDIAG(094) = -2378.812391459071;
HLSDIAG(095) = -2378.810768270247;
HLSDIAG(096) = -2378.807565118553;
HLSDIAG(097) = -2378.840045209658;
HLSDIAG(098) = -2378.837607801098;
HLSDIAG(099) = -2378.806900284714;
HLSDIAG(100) = -2378.815303666714;
HLSDIAG(101) = -2378.807757383343;
HLSDIAG(102) = -2378.819707746587;
HLSDIAG(103) = -2378.800854562318;
HLSDIAG(104) = -2378.793939517898;
HLSDIAG(105) = -2378.765365611082;
HLSDIAG(106) = -2378.835107851576;
HLSDIAG(107) = -2378.763994396071;
HLSDIAG(108) = -2378.833738437859;
HLSDIAG(109) = -2378.818285233128;
HLSDIAG(110) = -2378.816921806664;
HLSDIAG(111) = -2378.739433679493;
HLSDIAG(112) = -2378.739253602328;
HLSDIAG(113) = -2378.821308768619;
HLSDIAG(114) = -2378.821522623452;
HLSDIAG(115) = -2378.827224271811;
HLSDIAG(116) = -2378.826436697590;
HLSDIAG(117) = -2378.776782208817;
HLSDIAG(118) = -2378.774984858991;
HLSDIAG(119) = -2378.729712229449;
HLSDIAG(120) = -2378.729477413718;
HLSDIAG(121) = -2378.805783025024;
HLSDIAG(122) = -2378.804851031829;
 
 {CI;HLSMAT,LS,TRANLS,\
    6041.2,\               ! ^4B2g  wavefunction
    6051.2,\               ! ^4B3g  wavefunction
    6061.2,\               ! ^4B1u  wavefunction
    6071.2,                ! ^4Au   wavefunction
 
    print,\
    HLS = 1,\              ! The property matrices are printed
    VLS = 1;               ! The eigenvectors are printed

    option,\
    HLSTRANS = 0,\         ! SO matrix calculation should be performed in the not spin-symmetry adapted basis.
    MATEL = 1;}            ! The entire SO matrix is calculated using HLSMAT and individual matrix elements are printed







