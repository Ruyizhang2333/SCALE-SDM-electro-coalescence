	  �  =   k820309              15.0        ��`                                                                                                           
       admin/mod_init_driver.f90 MOD_INIT_DRIVER              SCALELES_INIT                      @                             
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
                                 &           #         @                                  '                     #         @                                   (                   #SCALELES_INIT%GAS_CTG )   #SCALELES_INIT%LKMAX *   #SCALELES_INIT%JMAX +   #SCALELES_INIT%IMAX ,   #SCALELES_INIT%QA -   #SCALELES_INIT%KMAX .   #SCALELES_INIT%JA /   #SCALELES_INIT%IA 0   #SCALELES_INIT%KA 1   #SCALELES_INIT%DAUGHTER_JA 2   #SCALELES_INIT%DAUGHTER_IA 3   #SCALELES_INIT%DAUGHTER_KA 4   #SCALELES_INIT%PARENT_JA 5   #SCALELES_INIT%PARENT_IA 6   #SCALELES_INIT%PARENT_KA 7   #COMM_WORLD 8   #INTERCOMM_PARENT 9   #INTERCOMM_CHILD :   #CNF_FNAME ;                                             )                     D                                *                     D                                +                     D                                ,                                                      -                     D                                .                                                      /                                                      0                                                      1                                                      2                         p          p            p                                                                    3                         p          p            p                                                                    4                         p          p            p                                                                    5                         p          p            p                                                                    6                         p          p            p                                                                    7                         p          p            p                                    
  @                               8                     
  @                               9                     
  @                               :                     
  @                              ;                              �   2      fn#fn %   �      b   uapp(MOD_INIT_DRIVER     �   @   J  SCALE_PRECISION    0  @   J  SCALE_STDIO    p  @   J  SCALE_PROF    �  H   J  DC_LOG    �  M   J  GTOOL_FILE    E  �       LOGINIT+DC_LOG /   �  @      LOGINIT%PRESENT+DC_LOG=PRESENT )     =      LOGINIT%TRIM+DC_LOG=TRIM (   X  @   a   LOGINIT%FID_CONF+DC_LOG '   �  @   a   LOGINIT%FID_LOG+DC_LOG &   �  @   a   LOGINIT%MASTER+DC_LOG (     H       FILECLOSEALL+GTOOL_FILE "   `  r       H_MID+SCALE_STDIO &   �  H       PROF_SETUP+SCALE_PROF #     s       H_LONG+SCALE_STDIO %   �  �       IO_SETUP+SCALE_STDIO 5   3  @      IO_SETUP%PRESENT+SCALE_STDIO=PRESENT /   s  =      IO_SETUP%TRIM+SCALE_STDIO=TRIM /   �  P   a   IO_SETUP%MODELNAME+SCALE_STDIO 8      @   a   IO_SETUP%CALL_FROM_LAUNCHER+SCALE_STDIO .   @  P   a   IO_SETUP%FNAME_IN+SCALE_STDIO )   �  z       IO_LOG_SETUP+SCALE_STDIO 3   
  =      IO_LOG_SETUP%TRIM+SCALE_STDIO=TRIM 0   G  @   a   IO_LOG_SETUP%MYRANK+SCALE_STDIO 3   �  @   a   IO_LOG_SETUP%IS_MASTER+SCALE_STDIO (   �  @       IO_FID_CONF+SCALE_STDIO '   	  @       IO_FID_LOG+SCALE_STDIO !   G	  @       IO_L+SCALE_STDIO )   �	  �       PROF_RAPSTART+SCALE_PROF 3   
  =      PROF_RAPSTART%TRIM+SCALE_PROF=TRIM 9   \
  @      PROF_RAPSTART%PRESENT+SCALE_PROF=PRESENT 6   �
  L   a   PROF_RAPSTART%RAPNAME_BASE+SCALE_PROF /   �
  @   a   PROF_RAPSTART%LEVEL+SCALE_PROF '   (  �       PROF_RAPEND+SCALE_PROF 1   �  =      PROF_RAPEND%TRIM+SCALE_PROF=TRIM 7   �  @      PROF_RAPEND%PRESENT+SCALE_PROF=PRESENT 4   9  L   a   PROF_RAPEND%RAPNAME_BASE+SCALE_PROF -   �  @   a   PROF_RAPEND%LEVEL+SCALE_PROF *   �  H       PROF_RAPREPORT+SCALE_PROF            SCALELES_INIT ;   '  @     SCALELES_INIT%GAS_CTG+SCALE_TRACER=GAS_CTG @   g  @     SCALELES_INIT%LKMAX+SCALE_LAND_GRID_INDEX=LKMAX 9   �  @     SCALELES_INIT%JMAX+SCALE_GRID_INDEX=JMAX 9   �  @     SCALELES_INIT%IMAX+SCALE_GRID_INDEX=IMAX 1   '  @     SCALELES_INIT%QA+SCALE_TRACER=QA 9   g  @     SCALELES_INIT%KMAX+SCALE_GRID_INDEX=KMAX 5   �  @     SCALELES_INIT%JA+SCALE_GRID_INDEX=JA 5   �  @     SCALELES_INIT%IA+SCALE_GRID_INDEX=IA 5   '  @     SCALELES_INIT%KA+SCALE_GRID_INDEX=KA F   g  �     SCALELES_INIT%DAUGHTER_JA+SCALE_GRID_NEST=DAUGHTER_JA F   �  �     SCALELES_INIT%DAUGHTER_IA+SCALE_GRID_NEST=DAUGHTER_IA F   �  �     SCALELES_INIT%DAUGHTER_KA+SCALE_GRID_NEST=DAUGHTER_KA B   #  �     SCALELES_INIT%PARENT_JA+SCALE_GRID_NEST=PARENT_JA B   �  �     SCALELES_INIT%PARENT_IA+SCALE_GRID_NEST=PARENT_IA B   K  �     SCALELES_INIT%PARENT_KA+SCALE_GRID_NEST=PARENT_KA )   �  @   a   SCALELES_INIT%COMM_WORLD /     @   a   SCALELES_INIT%INTERCOMM_PARENT .   _  @   a   SCALELES_INIT%INTERCOMM_CHILD (   �  P   a   SCALELES_INIT%CNF_FNAME 