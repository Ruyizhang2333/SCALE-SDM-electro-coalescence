	  aB  d   k820309              15.0        ��`                                                                                                           
       atmos-les/dynamics/scale_atmos_dyn_rk_heve.F90 SCALE_ATMOS_DYN_RK_HEVE              ATMOS_DYN_RK_HEVE_SETUP ATMOS_DYN_RK_HEVE                      @                             
                            @                             
                            @                              
                            @                              
                            @                             
                            @                              
                                                                                                                                                                    	                                                         
                                                                                                                                                                                                                                                                                                                                                             D                                                                                                                 D                                                                                                                                                                                                                                                                   1                                                                                                   1             @                                                                                                                                              2             @                                                                                                                                              3                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5#         @                                                       #ATMOS_DYN_TYPE     #BND_W !   #BND_E "   #BND_S #   #BND_N $             
                                                      1           
                                  !                     
                                  "                     
                                  #                     
                                  $           #         @                                   %               7   #ATMOS_DYN_RK_HEVE%JA &   #ATMOS_DYN_RK_HEVE%KA '   #ATMOS_DYN_RK_HEVE%IA (   #ATMOS_DYN_RK_HEVE%SIGN )   #ATMOS_DYN_RK_HEVE%ABS *   #ATMOS_DYN_RK_HEVE%MIN +   #DENS_RK ,   #MOMZ_RK -   #MOMX_RK .   #MOMY_RK /   #RHOT_RK 0   #MFLX_HI 1   #TFLX_HI 2   #DENS0 3   #MOMZ0 4   #MOMX0 5   #MOMY0 6   #RHOT0 7   #DENS 8   #MOMZ 9   #MOMX :   #MOMY ;   #RHOT <   #DENS_T =   #MOMZ_T >   #MOMX_T ?   #MOMY_T @   #RHOT_T A   #RTOT B   #CVTOT C   #CORIOLI D   #NUM_DIFF E   #DIVDMP_COEF F   #FLAG_FCT_RHO G   #FLAG_FCT_MOMENTUM H   #FLAG_FCT_T I   #FLAG_FCT_ALONG_STREAM J   #CDZ K   #FDZ L   #FDX M   #FDY N   #RCDZ O   #RCDX P   #RCDY Q   #RFDZ R   #RFDX S   #RFDY T   #PHI U   #GSQRT V   #J13G W   #J23G X   #J33G Y   #MAPF Z   #REF_PRES [   #REF_DENS \   #BND_W ]   #BND_E ^   #BND_S _   #BND_N `   #DTRK a   #DT b                                                                                                         &                                                      '                                                      (                          @                           )     SIGN               @                           *     ABS               @                           +     MIN          D @                              ,                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                -                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                .                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                /                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               D                                0                    
         p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
D @                              1                    
           p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                   D @                              2                    
           p        5 r    p        5 r    p        5 r    p          5 r      5 r      5 r      p            5 r      5 r      5 r      p                                   
  `                              3                    
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
  @                               J                    
                                 K                    
    p          5 r        5 r                               
                                 L                    
    p           5 r    n                                       1     5 r    n                                      1                                    
                                 M                    
    p           5 r    n                                       1     5 r    n                                      1                                    
                                 N                    
    p           5 r    n                                       1     5 r    n                                      1                                    
  @                              O                    
    p          5 r        5 r                               
  @                              P                    
     p          5 r        5 r                               
  @                              Q                    
 !   p          5 r        5 r                               
  @                              R                    
 "   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              S                    
 #   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              T                    
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
  @                              a     
                
                                 b     
         �   O      fn#fn -   �   :   b   uapp(SCALE_ATMOS_DYN_RK_HEVE     )  @   J  SCALE_PRECISION    i  @   J  SCALE_STDIO    �  @   J  SCALE_PROF !   �  @   j  SCALE_GRID_INDEX    )  @   J  SCALE_INDEX    i  @   J  SCALE_TRACER $   �  @       IS+SCALE_GRID_INDEX !   �  @       IO_L+SCALE_STDIO '   )  @       IO_FID_LOG+SCALE_STDIO #   i  p       RP+SCALE_PRECISION $   �  @       KA+SCALE_GRID_INDEX $     @       IA+SCALE_GRID_INDEX $   Y  @       JA+SCALE_GRID_INDEX $   �  @       JS+SCALE_GRID_INDEX $   �  @       JE+SCALE_GRID_INDEX (     @       JBLOCK+SCALE_GRID_INDEX $   Y  @       IE+SCALE_GRID_INDEX (   �  @       IBLOCK+SCALE_GRID_INDEX $   �  @       KS+SCALE_GRID_INDEX $     @       KE+SCALE_GRID_INDEX &   Y  q       ZDIR+SCALE_GRID_INDEX #   �  q       I_DENS+SCALE_INDEX %   ;  @       IEH+SCALE_GRID_INDEX &   {  q       XDIR+SCALE_GRID_INDEX %   �  @       JEH+SCALE_GRID_INDEX &   ,  q       YDIR+SCALE_GRID_INDEX #   �  q       I_MOMZ+SCALE_INDEX #   	  q       I_MOMX+SCALE_INDEX #   	  q       I_MOMY+SCALE_INDEX #   �	  q       I_RHOT+SCALE_INDEX (   a
  �       ATMOS_DYN_RK_HEVE_SETUP 7   �
  L   a   ATMOS_DYN_RK_HEVE_SETUP%ATMOS_DYN_TYPE .   5  @   a   ATMOS_DYN_RK_HEVE_SETUP%BND_W .   u  @   a   ATMOS_DYN_RK_HEVE_SETUP%BND_E .   �  @   a   ATMOS_DYN_RK_HEVE_SETUP%BND_S .   �  @   a   ATMOS_DYN_RK_HEVE_SETUP%BND_N "   5  �      ATMOS_DYN_RK_HEVE 9   �  @     ATMOS_DYN_RK_HEVE%JA+SCALE_GRID_INDEX=JA 9   $  @     ATMOS_DYN_RK_HEVE%KA+SCALE_GRID_INDEX=KA 9   d  @     ATMOS_DYN_RK_HEVE%IA+SCALE_GRID_INDEX=IA '   �  =      ATMOS_DYN_RK_HEVE%SIGN &   �  <      ATMOS_DYN_RK_HEVE%ABS &     <      ATMOS_DYN_RK_HEVE%MIN *   Y    a   ATMOS_DYN_RK_HEVE%DENS_RK *   m    a   ATMOS_DYN_RK_HEVE%MOMZ_RK *   �    a   ATMOS_DYN_RK_HEVE%MOMX_RK *   �    a   ATMOS_DYN_RK_HEVE%MOMY_RK *   �    a   ATMOS_DYN_RK_HEVE%RHOT_RK *   �  T  a   ATMOS_DYN_RK_HEVE%MFLX_HI *     T  a   ATMOS_DYN_RK_HEVE%TFLX_HI (   e    a   ATMOS_DYN_RK_HEVE%DENS0 (   y    a   ATMOS_DYN_RK_HEVE%MOMZ0 (   �    a   ATMOS_DYN_RK_HEVE%MOMX0 (   �    a   ATMOS_DYN_RK_HEVE%MOMY0 (   �    a   ATMOS_DYN_RK_HEVE%RHOT0 '   �    a   ATMOS_DYN_RK_HEVE%DENS '   �    a   ATMOS_DYN_RK_HEVE%MOMZ '   �     a   ATMOS_DYN_RK_HEVE%MOMX '   "    a   ATMOS_DYN_RK_HEVE%MOMY '   #    a   ATMOS_DYN_RK_HEVE%RHOT )   -$    a   ATMOS_DYN_RK_HEVE%DENS_T )   A%    a   ATMOS_DYN_RK_HEVE%MOMZ_T )   U&    a   ATMOS_DYN_RK_HEVE%MOMX_T )   i'    a   ATMOS_DYN_RK_HEVE%MOMY_T )   }(    a   ATMOS_DYN_RK_HEVE%RHOT_T '   �)    a   ATMOS_DYN_RK_HEVE%RTOT (   �*    a   ATMOS_DYN_RK_HEVE%CVTOT *   �+    a   ATMOS_DYN_RK_HEVE%CORIOLI +   �,  �  a   ATMOS_DYN_RK_HEVE%NUM_DIFF .   a.  @   a   ATMOS_DYN_RK_HEVE%DIVDMP_COEF /   �.  @   a   ATMOS_DYN_RK_HEVE%FLAG_FCT_RHO 4   �.  @   a   ATMOS_DYN_RK_HEVE%FLAG_FCT_MOMENTUM -   !/  @   a   ATMOS_DYN_RK_HEVE%FLAG_FCT_T 8   a/  @   a   ATMOS_DYN_RK_HEVE%FLAG_FCT_ALONG_STREAM &   �/  �   a   ATMOS_DYN_RK_HEVE%CDZ &   50    a   ATMOS_DYN_RK_HEVE%FDZ &   ;1    a   ATMOS_DYN_RK_HEVE%FDX &   A2    a   ATMOS_DYN_RK_HEVE%FDY '   G3  �   a   ATMOS_DYN_RK_HEVE%RCDZ '   �3  �   a   ATMOS_DYN_RK_HEVE%RCDX '   o4  �   a   ATMOS_DYN_RK_HEVE%RCDY '   5    a   ATMOS_DYN_RK_HEVE%RFDZ '   	6    a   ATMOS_DYN_RK_HEVE%RFDX '   7    a   ATMOS_DYN_RK_HEVE%RFDY &   8    a   ATMOS_DYN_RK_HEVE%PHI (   )9  T  a   ATMOS_DYN_RK_HEVE%GSQRT '   }:  T  a   ATMOS_DYN_RK_HEVE%J13G '   �;  T  a   ATMOS_DYN_RK_HEVE%J23G '   %=  @   a   ATMOS_DYN_RK_HEVE%J33G '   e=  T  a   ATMOS_DYN_RK_HEVE%MAPF +   �>    a   ATMOS_DYN_RK_HEVE%REF_PRES +   �?    a   ATMOS_DYN_RK_HEVE%REF_DENS (   �@  @   a   ATMOS_DYN_RK_HEVE%BND_W (   !A  @   a   ATMOS_DYN_RK_HEVE%BND_E (   aA  @   a   ATMOS_DYN_RK_HEVE%BND_S (   �A  @   a   ATMOS_DYN_RK_HEVE%BND_N '   �A  @   a   ATMOS_DYN_RK_HEVE%DTRK %   !B  @   a   ATMOS_DYN_RK_HEVE%DT 