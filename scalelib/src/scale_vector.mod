	  2   M   k820309              15.0        ��`                                                                                                           
       common/scale_vector.F90 SCALE_VECTOR       
       VECTR_XYZ2LATLON VECTR_LATLON2XYZ VECTR_CROSS VECTR_DOT VECTR_ABS VECTR_ANGLE VECTR_INTERSEC VECTR_ANTICLOCKWISE VECTR_TRIANGLE VECTR_DISTANCE                      @                             
                            @                             
                            @                             
                                                                                                             #         @                                                      #VECTR_XYZ2LATLON%ACOS    #VECTR_XYZ2LATLON%ASIN    #VECTR_XYZ2LATLON%SQRT    #X 	   #Y 
   #Z    #LAT    #LON                  @                                ACOS               @                                ASIN               @                                SQRT           
                                 	     
                
                                 
     
                
                                      
                D                                     
                 D                                     
       #         @                                                      #VECTR_LATLON2XYZ%SIN    #VECTR_LATLON2XYZ%COS    #LAT    #LON    #X    #Y    #Z    #RADIUS                  @                                SIN               @                                COS           
  @                                   
                
  @                                   
                D                                     
                 D                                     
                 D                                     
                 
                                      
      #         @                                                      #NV    #A    #B    #C    #D              D                                                   
     p          p            p                                    
                                                    
    p          p            p                                    
                                                    
    p          p            p                                    
                                                    
    p          p            p                                    
                                                    
    p          p            p                          #         @                                                      #L    #A    #B     #C !   #D "             D                                     
                 
                                                    
    p          p            p                                    
                                                     
    p          p            p                                    
                                 !                   
    p          p            p                                    
                                 "                   
 	   p          p            p                          #         @                                  #                   #VECTR_ABS%SQRT $   #L %   #A &                 @                           $     SQRT           D @                              %     
                 
                                 &                   
 
   p          p            p                          #         @                                  '                   #VECTR_ANGLE%ATAN2 (   #ANGLE )   #A *   #B +   #C ,                 @                           (     ATAN2           D                                )     
                 
  @                              *                   
    p          p            p                                    
  @                              +                   
    p          p            p                                    
  @                              ,                   
    p          p            p                          #         @                                   -                   #VECTR_INTERSEC%ABS .   #VECTR_INTERSEC%SIGN /   #IFCROSS 0   #P 1   #A 2   #B 3   #C 4   #D 5                 @                           .     ABS               @                           /     SIGN           D                                 0                      D @                              1                   
     p          p            p                                    
  @                              2                   
    p          p            p                                    
  @                              3                   
    p          p            p                                    
  @                              4                   
    p          p            p                                    
  @                              5                   
    p          p            p                          #         @                                   6                   #VECTR_ANTICLOCKWISE%ABS 7   #VERTEX 8   #NVERT 9                 @                           7     ABS          
D                                8                    
       p        5 � p        r 9   p          5 � p        r 9     p            5 � p        r 9     p                                    
                                  9           %         @                                 :                   
       #VECTR_TRIANGLE%O ;   #A <   #B =   #C >   #POLYGON_TYPE ?   #RADIUS @                 @                             ;                   
                                                    T
W
p          n
      
                                 0.0  h   p          p          p            p                                                        
                                 <                   
    p          p            p                                    
                                 =                   
    p          p            p                                    
                                 >                   
     p          p            p                                    
                                 ?                    1           
                                 @     
      #         @                                   A                   #VECTR_DISTANCE%ATAN2 B   #VECTR_DISTANCE%SIN C   #VECTR_DISTANCE%COS D   #VECTR_DISTANCE%SQRT E   #R F   #LON1 G   #LAT1 H   #LON2 I   #LAT2 J   #DIST K                 @                           B     ATAN2               @                           C     SIN               @                           D     COS               @                           E     SQRT           
                                 F     
                
                                 G     
                
  @                              H     
                
                                 I     
                
  @                              J     
                D                                K     
          �   -      fn#fn "   �   �   b   uapp(SCALE_VECTOR     l  @   J  SCALE_PRECISION    �  @   J  SCALE_STDIO    �  @   J  SCALE_PROF #   ,  p       RP+SCALE_PRECISION !   �  �       VECTR_XYZ2LATLON &   \  =      VECTR_XYZ2LATLON%ACOS &   �  =      VECTR_XYZ2LATLON%ASIN &   �  =      VECTR_XYZ2LATLON%SQRT #     @   a   VECTR_XYZ2LATLON%X #   S  @   a   VECTR_XYZ2LATLON%Y #   �  @   a   VECTR_XYZ2LATLON%Z %   �  @   a   VECTR_XYZ2LATLON%LAT %     @   a   VECTR_XYZ2LATLON%LON !   S  �       VECTR_LATLON2XYZ %     <      VECTR_LATLON2XYZ%SIN %   >  <      VECTR_LATLON2XYZ%COS %   z  @   a   VECTR_LATLON2XYZ%LAT %   �  @   a   VECTR_LATLON2XYZ%LON #   �  @   a   VECTR_LATLON2XYZ%X #   :  @   a   VECTR_LATLON2XYZ%Y #   z  @   a   VECTR_LATLON2XYZ%Z (   �  @   a   VECTR_LATLON2XYZ%RADIUS    �  l       VECTR_CROSS    f  �   a   VECTR_CROSS%NV    �  �   a   VECTR_CROSS%A    �	  �   a   VECTR_CROSS%B    "
  �   a   VECTR_CROSS%C    �
  �   a   VECTR_CROSS%D    J  k       VECTR_DOT    �  @   a   VECTR_DOT%L    �  �   a   VECTR_DOT%A    �  �   a   VECTR_DOT%B      �   a   VECTR_DOT%C    �  �   a   VECTR_DOT%D    E  j       VECTR_ABS    �  =      VECTR_ABS%SQRT    �  @   a   VECTR_ABS%L    ,  �   a   VECTR_ABS%A    �         VECTR_ANGLE "   ?  >      VECTR_ANGLE%ATAN2 "   }  @   a   VECTR_ANGLE%ANGLE    �  �   a   VECTR_ANGLE%A    Q  �   a   VECTR_ANGLE%B    �  �   a   VECTR_ANGLE%C    y  �       VECTR_INTERSEC #   "  <      VECTR_INTERSEC%ABS $   ^  =      VECTR_INTERSEC%SIGN '   �  @   a   VECTR_INTERSEC%IFCROSS !   �  �   a   VECTR_INTERSEC%P !   o  �   a   VECTR_INTERSEC%A !     �   a   VECTR_INTERSEC%B !   �  �   a   VECTR_INTERSEC%C !   +  �   a   VECTR_INTERSEC%D $   �  |       VECTR_ANTICLOCKWISE (   ;  <      VECTR_ANTICLOCKWISE%ABS +   w    a   VECTR_ANTICLOCKWISE%VERTEX *   {  @   a   VECTR_ANTICLOCKWISE%NVERT    �  �       VECTR_TRIANGLE !   T  ?     VECTR_TRIANGLE%O !   �  �   a   VECTR_TRIANGLE%A !   '  �   a   VECTR_TRIANGLE%B !   �  �   a   VECTR_TRIANGLE%C ,   O  L   a   VECTR_TRIANGLE%POLYGON_TYPE &   �  @   a   VECTR_TRIANGLE%RADIUS    �  �       VECTR_DISTANCE %   �  >      VECTR_DISTANCE%ATAN2 #   �  <      VECTR_DISTANCE%SIN #   9  <      VECTR_DISTANCE%COS $   u  =      VECTR_DISTANCE%SQRT !   �  @   a   VECTR_DISTANCE%R $   �  @   a   VECTR_DISTANCE%LON1 $   2  @   a   VECTR_DISTANCE%LAT1 $   r  @   a   VECTR_DISTANCE%LON2 $   �  @   a   VECTR_DISTANCE%LAT2 $   �  @   a   VECTR_DISTANCE%DIST 