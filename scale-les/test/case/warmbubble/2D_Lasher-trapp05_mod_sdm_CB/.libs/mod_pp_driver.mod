	  r  <   k820309              15.0        j$_                                                                                                           
       admin/mod_pp_driver.f90 MOD_PP_DRIVER              SCALELES_PP                      @                             
                            @                              
                            @                              
                            @                              
       LOGINIT                      @                              
       FILECLOSEALL #         @                                                     #LOGINIT%PRESENT    #LOGINIT%TRIM    #FID_CONF 	   #FID_LOG 
   #MASTER                  @                                PRESENT               @                                TRIM           
                                  	                     
                                 
                     
                                            #         @                                                                                                                                            @               64#         @                                                                                                                                                            256#         @                                                     #IO_SETUP%PRESENT    #IO_SETUP%TRIM    #MODELNAME    #CALL_FROM_LAUNCHER    #FNAME_IN                  @                                PRESENT               @                                TRIM           
                                      @                                
                                                       
                                                           #         @                                                     #IO_LOG_SETUP%TRIM    #MYRANK    #IS_MASTER                  @                                TRIM           
                                                       
                                                         @                                                       @                                                       @                                           #         @                                                     #PROF_RAPSTART%TRIM    #PROF_RAPSTART%PRESENT    #RAPNAME_BASE     #LEVEL !                 @                                TRIM               @                                PRESENT           
                                                      1           
                                 !           #         @                                  "                   #PROF_RAPEND%TRIM #   #PROF_RAPEND%PRESENT $   #RAPNAME_BASE %   #LEVEL &                 @                           #     TRIM               @                           $     PRESENT           
                                 %                    1           
                                 &           #         @                                  '                     #         @                                   (                   #SCALELES_PP%LKMAX )   #SCALELES_PP%JMAX *   #SCALELES_PP%IMAX +   #SCALELES_PP%QA ,   #SCALELES_PP%KMAX -   #SCALELES_PP%JA .   #SCALELES_PP%IA /   #SCALELES_PP%KA 0   #SCALELES_PP%DAUGHTER_JA 1   #SCALELES_PP%DAUGHTER_IA 2   #SCALELES_PP%DAUGHTER_KA 3   #SCALELES_PP%PARENT_JA 4   #SCALELES_PP%PARENT_IA 5   #SCALELES_PP%PARENT_KA 6   #COMM_WORLD 7   #INTERCOMM_PARENT 8   #INTERCOMM_CHILD 9   #CNF_FNAME :            D                                )                     D                                *                     D                                +                                                      ,                     D                                -                                                      .                                                      /                                                      0                                                      1                         p          p            p                                                                    2                         p          p            p                                                                    3                         p          p            p                                                                    4                         p          p            p                                                                    5                         p          p            p                                                                    6                         p          p            p                                    
  @                               7                     
  @                               8                     
  @                               9                     
  @                              :                              �   .      fn#fn #   �      b   uapp(MOD_PP_DRIVER     �   @   J  SCALE_PRECISION    *  @   J  SCALE_STDIO    j  @   J  SCALE_PROF    �  H   J  DC_LOG    �  M   J  GTOOL_FILE    ?  �       LOGINIT+DC_LOG /   �  @      LOGINIT%PRESENT+DC_LOG=PRESENT )     =      LOGINIT%TRIM+DC_LOG=TRIM (   R  @   a   LOGINIT%FID_CONF+DC_LOG '   �  @   a   LOGINIT%FID_LOG+DC_LOG &   �  @   a   LOGINIT%MASTER+DC_LOG (     H       FILECLOSEALL+GTOOL_FILE "   Z  r       H_MID+SCALE_STDIO &   �  H       PROF_SETUP+SCALE_PROF #     s       H_LONG+SCALE_STDIO %   �  �       IO_SETUP+SCALE_STDIO 5   -  @      IO_SETUP%PRESENT+SCALE_STDIO=PRESENT /   m  =      IO_SETUP%TRIM+SCALE_STDIO=TRIM /   �  P   a   IO_SETUP%MODELNAME+SCALE_STDIO 8   �  @   a   IO_SETUP%CALL_FROM_LAUNCHER+SCALE_STDIO .   :  P   a   IO_SETUP%FNAME_IN+SCALE_STDIO )   �  z       IO_LOG_SETUP+SCALE_STDIO 3     =      IO_LOG_SETUP%TRIM+SCALE_STDIO=TRIM 0   A  @   a   IO_LOG_SETUP%MYRANK+SCALE_STDIO 3   �  @   a   IO_LOG_SETUP%IS_MASTER+SCALE_STDIO (   �  @       IO_FID_CONF+SCALE_STDIO '   	  @       IO_FID_LOG+SCALE_STDIO !   A	  @       IO_L+SCALE_STDIO )   �	  �       PROF_RAPSTART+SCALE_PROF 3   
  =      PROF_RAPSTART%TRIM+SCALE_PROF=TRIM 9   V
  @      PROF_RAPSTART%PRESENT+SCALE_PROF=PRESENT 6   �
  L   a   PROF_RAPSTART%RAPNAME_BASE+SCALE_PROF /   �
  @   a   PROF_RAPSTART%LEVEL+SCALE_PROF '   "  �       PROF_RAPEND+SCALE_PROF 1   �  =      PROF_RAPEND%TRIM+SCALE_PROF=TRIM 7   �  @      PROF_RAPEND%PRESENT+SCALE_PROF=PRESENT 4   3  L   a   PROF_RAPEND%RAPNAME_BASE+SCALE_PROF -     @   a   PROF_RAPEND%LEVEL+SCALE_PROF *   �  H       PROF_RAPREPORT+SCALE_PROF      �      SCALELES_PP >   �  @     SCALELES_PP%LKMAX+SCALE_LAND_GRID_INDEX=LKMAX 7   *  @     SCALELES_PP%JMAX+SCALE_GRID_INDEX=JMAX 7   j  @     SCALELES_PP%IMAX+SCALE_GRID_INDEX=IMAX /   �  @     SCALELES_PP%QA+SCALE_TRACER=QA 7   �  @     SCALELES_PP%KMAX+SCALE_GRID_INDEX=KMAX 3   *  @     SCALELES_PP%JA+SCALE_GRID_INDEX=JA 3   j  @     SCALELES_PP%IA+SCALE_GRID_INDEX=IA 3   �  @     SCALELES_PP%KA+SCALE_GRID_INDEX=KA D   �  �     SCALELES_PP%DAUGHTER_JA+SCALE_GRID_NEST=DAUGHTER_JA D   ~  �     SCALELES_PP%DAUGHTER_IA+SCALE_GRID_NEST=DAUGHTER_IA D     �     SCALELES_PP%DAUGHTER_KA+SCALE_GRID_NEST=DAUGHTER_KA @   �  �     SCALELES_PP%PARENT_JA+SCALE_GRID_NEST=PARENT_JA @   :  �     SCALELES_PP%PARENT_IA+SCALE_GRID_NEST=PARENT_IA @   �  �     SCALELES_PP%PARENT_KA+SCALE_GRID_NEST=PARENT_KA '   b  @   a   SCALELES_PP%COMM_WORLD -   �  @   a   SCALELES_PP%INTERCOMM_PARENT ,   �  @   a   SCALELES_PP%INTERCOMM_CHILD &   "  P   a   SCALELES_PP%CNF_FNAME 