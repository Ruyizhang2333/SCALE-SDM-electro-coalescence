	  fB  d   k820309              15.0        ��`                                                                                                           
       atmos-les/dynamics/scale_atmos_dyn_rk_hevi.F90 SCALE_ATMOS_DYN_RK_HEVI              ATMOS_DYN_RK_HEVI_SETUP ATMOS_DYN_RK_HEVI                      @                             
                            @                             
                            @                              
                            @                              
                            @                             
                            @                              
                                                                                                                                                                                                                  	                                                        
                                                                                                                                                                                             D                                                                                                                                                                          D                                                                                                                 D                                                                                                                                                                                                                                                                   1                                                                                                   1             @                                                                                                                                              2             @                                                                                                                                              3                                                                                                   2                                                                                                   5                                                                                                   3                                                                                                   4#         @                                                        #ATMOS_DYN_TYPE !   #BND_W "   #BND_E #   #BND_S $   #BND_N %             
                                 !                    1           
                                  "                     
                                  #                     
                                  $                     
                                  %           #         @                                   &               7   #ATMOS_DYN_RK_HEVI%JA '   #ATMOS_DYN_RK_HEVI%IA (   #ATMOS_DYN_RK_HEVI%KA )   #ATMOS_DYN_RK_HEVI%MIN *   #ATMOS_DYN_RK_HEVI%MAX +   #DENS_RK ,   #MOMZ_RK -   #MOMX_RK .   #MOMY_RK /   #RHOT_RK 0   #MFLX_HI 1   #TFLX_HI 2   #DENS0 3   #MOMZ0 4   #MOMX0 5   #MOMY0 6   #RHOT0 7   #DENS 8   #MOMZ 9   #MOMX :   #MOMY ;   #RHOT <   #DENS_T =   #MOMZ_T >   #MOMX_T ?   #MOMY_T @   #RHOT_T A   #RTOT B   #CVTOT C   #CORIOLI D   #NUM_DIFF E   #DIVDMP_COEF F   #FLAG_FCT_RHO G   #FLAG_FCT_MOMENTUM H   #FLAG_FCT_T I   #FLAG_FCT_ALONG_STREAM J   #CDZ K   #FDZ L   #FDX M   #FDY N   #RCDZ O   #RCDX P   #RCDY Q   #RFDZ R   #RFDX S   #RFDY T   #PHI U   #GSQRT V   #J13G W   #J23G X   #J33G Y   #MAPF Z   #REF_PRES [   #REF_DENS \   #BND_W ]   #BND_E ^   #BND_S _   #BND_N `   #DTRK a   #DT b                                                                                                                                       '                                                      (                                                      )                          @                           *     MIN               @                           +     MAX          D                                ,                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                -                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                .                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                /                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                0                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
D                                1                    
           p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                   D                                2                    
           p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                   
                                 3                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 4                    
 	       p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 5                    
 
       p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 6                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 7                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 8                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 9                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 :                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 ;                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 <                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 =                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 >                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 ?                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 @                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 A                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 B                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 C                    
        p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 D                    
        p        5 r    p        p        p          p          5 r      5 r        p          5 r      5 r                               
                                 E                    
            p        p        p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p          p            5 r      5 r      5 r      p          p                                    
                                 F     
                
                                  G                     
                                  H                     
                                  I                     
                                  J                    
                                 K                    
    p          5 r        5 r                               
                                 L                    
    p           5 r    n                                       1     5 r    n                                      1                                    
                                 M                    
    p           5 r    n                                       1     5 r    n                                      1                                    
                                 N                    
    p           5 r    n                                       1     5 r    n                                      1                                    
                                 O                    
    p          5 r        5 r                               
                                 P                    
     p          5 r        5 r                               
                                 Q                    
 !   p          5 r        5 r                               
                                 R                    
 "   p           5 r    n                                       1     5 r    n                                      1                                    
                                 S                    
 #   p           5 r    n                                       1     5 r    n                                      1                                    
                                 T                    
 $   p           5 r    n                                       1     5 r    n                                      1                                    
                                 U                    
 %       p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 V                    
 &         p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                   
                                 W                    
 '         p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                   
                                 X                    
 (         p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                    
                                 Y     
               
                                 Z                    
 )         p        p        p        5 r    p        5 r    p          5 r      5 r      p          p            5 r      5 r      p          p                                   
                                 [                    
 *       p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                 \                    
 +       p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                                
                                  ]                     
                                  ^                     
                                  _                     
                                  `                     
                                 a     
                
                                 b     
         �   O      fn#fn -   �   :   b   uapp(SCALE_ATMOS_DYN_RK_HEVI     )  @   J  SCALE_PRECISION    i  @   J  SCALE_STDIO    �  @   J  SCALE_PROF !   �  @   j  SCALE_GRID_INDEX    )  @   J  SCALE_INDEX    i  @   J  SCALE_TRACER $   �  @       IS+SCALE_GRID_INDEX #   �  p       RP+SCALE_PRECISION !   Y  @       IO_L+SCALE_STDIO '   �  @       IO_FID_LOG+SCALE_STDIO $   �  @       KA+SCALE_GRID_INDEX $     @       IA+SCALE_GRID_INDEX $   Y  @       JA+SCALE_GRID_INDEX &   �  @       KMAX+SCALE_GRID_INDEX $   �  @       JS+SCALE_GRID_INDEX $     @       JE+SCALE_GRID_INDEX (   Y  @       JBLOCK+SCALE_GRID_INDEX $   �  @       IE+SCALE_GRID_INDEX (   �  @       IBLOCK+SCALE_GRID_INDEX $     @       KS+SCALE_GRID_INDEX $   Y  @       KE+SCALE_GRID_INDEX &   �  q       ZDIR+SCALE_GRID_INDEX #   
  q       I_DENS+SCALE_INDEX %   {  @       IEH+SCALE_GRID_INDEX &   �  q       XDIR+SCALE_GRID_INDEX %   ,  @       JEH+SCALE_GRID_INDEX &   l  q       YDIR+SCALE_GRID_INDEX #   �  q       I_MOMZ+SCALE_INDEX #   N	  q       I_RHOT+SCALE_INDEX #   �	  q       I_MOMX+SCALE_INDEX #   0
  q       I_MOMY+SCALE_INDEX (   �
  �       ATMOS_DYN_RK_HEVI_SETUP 7   )  L   a   ATMOS_DYN_RK_HEVI_SETUP%ATMOS_DYN_TYPE .   u  @   a   ATMOS_DYN_RK_HEVI_SETUP%BND_W .   �  @   a   ATMOS_DYN_RK_HEVI_SETUP%BND_E .   �  @   a   ATMOS_DYN_RK_HEVI_SETUP%BND_S .   5  @   a   ATMOS_DYN_RK_HEVI_SETUP%BND_N "   u  �      ATMOS_DYN_RK_HEVI 9   &  @     ATMOS_DYN_RK_HEVI%JA+SCALE_GRID_INDEX=JA 9   f  @     ATMOS_DYN_RK_HEVI%IA+SCALE_GRID_INDEX=IA 9   �  @     ATMOS_DYN_RK_HEVI%KA+SCALE_GRID_INDEX=KA &   �  <      ATMOS_DYN_RK_HEVI%MIN &   "  <      ATMOS_DYN_RK_HEVI%MAX *   ^    a   ATMOS_DYN_RK_HEVI%DENS_RK *   r    a   ATMOS_DYN_RK_HEVI%MOMZ_RK *   �    a   ATMOS_DYN_RK_HEVI%MOMX_RK *   �    a   ATMOS_DYN_RK_HEVI%MOMY_RK *   �    a   ATMOS_DYN_RK_HEVI%RHOT_RK *   �  T  a   ATMOS_DYN_RK_HEVI%MFLX_HI *     T  a   ATMOS_DYN_RK_HEVI%TFLX_HI (   j    a   ATMOS_DYN_RK_HEVI%DENS0 (   ~    a   ATMOS_DYN_RK_HEVI%MOMZ0 (   �    a   ATMOS_DYN_RK_HEVI%MOMX0 (   �    a   ATMOS_DYN_RK_HEVI%MOMY0 (   �    a   ATMOS_DYN_RK_HEVI%RHOT0 '   �    a   ATMOS_DYN_RK_HEVI%DENS '   �    a   ATMOS_DYN_RK_HEVI%MOMZ '   �     a   ATMOS_DYN_RK_HEVI%MOMX '   
"    a   ATMOS_DYN_RK_HEVI%MOMY '   #    a   ATMOS_DYN_RK_HEVI%RHOT )   2$    a   ATMOS_DYN_RK_HEVI%DENS_T )   F%    a   ATMOS_DYN_RK_HEVI%MOMZ_T )   Z&    a   ATMOS_DYN_RK_HEVI%MOMX_T )   n'    a   ATMOS_DYN_RK_HEVI%MOMY_T )   �(    a   ATMOS_DYN_RK_HEVI%RHOT_T '   �)    a   ATMOS_DYN_RK_HEVI%RTOT (   �*    a   ATMOS_DYN_RK_HEVI%CVTOT *   �+    a   ATMOS_DYN_RK_HEVI%CORIOLI +   �,  �  a   ATMOS_DYN_RK_HEVI%NUM_DIFF .   f.  @   a   ATMOS_DYN_RK_HEVI%DIVDMP_COEF /   �.  @   a   ATMOS_DYN_RK_HEVI%FLAG_FCT_RHO 4   �.  @   a   ATMOS_DYN_RK_HEVI%FLAG_FCT_MOMENTUM -   &/  @   a   ATMOS_DYN_RK_HEVI%FLAG_FCT_T 8   f/  @   a   ATMOS_DYN_RK_HEVI%FLAG_FCT_ALONG_STREAM &   �/  �   a   ATMOS_DYN_RK_HEVI%CDZ &   :0    a   ATMOS_DYN_RK_HEVI%FDZ &   @1    a   ATMOS_DYN_RK_HEVI%FDX &   F2    a   ATMOS_DYN_RK_HEVI%FDY '   L3  �   a   ATMOS_DYN_RK_HEVI%RCDZ '   �3  �   a   ATMOS_DYN_RK_HEVI%RCDX '   t4  �   a   ATMOS_DYN_RK_HEVI%RCDY '   5    a   ATMOS_DYN_RK_HEVI%RFDZ '   6    a   ATMOS_DYN_RK_HEVI%RFDX '   7    a   ATMOS_DYN_RK_HEVI%RFDY &   8    a   ATMOS_DYN_RK_HEVI%PHI (   .9  T  a   ATMOS_DYN_RK_HEVI%GSQRT '   �:  T  a   ATMOS_DYN_RK_HEVI%J13G '   �;  T  a   ATMOS_DYN_RK_HEVI%J23G '   *=  @   a   ATMOS_DYN_RK_HEVI%J33G '   j=  T  a   ATMOS_DYN_RK_HEVI%MAPF +   �>    a   ATMOS_DYN_RK_HEVI%REF_PRES +   �?    a   ATMOS_DYN_RK_HEVI%REF_DENS (   �@  @   a   ATMOS_DYN_RK_HEVI%BND_W (   &A  @   a   ATMOS_DYN_RK_HEVI%BND_E (   fA  @   a   ATMOS_DYN_RK_HEVI%BND_S (   �A  @   a   ATMOS_DYN_RK_HEVI%BND_N '   �A  @   a   ATMOS_DYN_RK_HEVI%DTRK %   &B  @   a   ATMOS_DYN_RK_HEVI%DT 