	  S  ?   k820309              15.0        P_                                                                                                           
       admin/mod_les_driver.f90 MOD_LES_DRIVER              SCALELES                      @                             
                            @                              
                            @                              
                            @                              
       LOGINIT                      @                              
       FILECLOSEALL #         @                                                     #LOGINIT%PRESENT    #LOGINIT%TRIM    #FID_CONF 	   #FID_LOG 
   #MASTER                  @                                PRESENT               @                                TRIM           
                                  	                     
                                 
                     
                                            #         @                                                                                                                                            @               64                                                                                                    256#         @                                                     #IO_SETUP%PRESENT    #IO_SETUP%TRIM    #MODELNAME    #CALL_FROM_LAUNCHER    #FNAME_IN                  @                                PRESENT               @                                TRIM           
                                      @                                
                                                       
                                                           #         @                                                     #IO_LOG_SETUP%TRIM    #MYRANK    #IS_MASTER                  @                                TRIM           
                                                       
                                                         @                                                       @                                                       @                                           #         @                                                       #         @                                                     #PROF_SETPREFX%TRIM    #PREFXNAME                  @                                TRIM           
                                                     1 #         @                                                      #PROF_RAPSTART%TRIM !   #PROF_RAPSTART%PRESENT "   #RAPNAME_BASE #   #LEVEL $                 @                           !     TRIM               @                           "     PRESENT           
                                 #                    1           
                                 $           #         @                                  %                   #PROF_RAPEND%TRIM &   #PROF_RAPEND%PRESENT '   #RAPNAME_BASE (   #LEVEL )                 @                           &     TRIM               @                           '     PRESENT           
                                 (                    1           
                                 )           #         @                                  *                     #         @                                   +                   #SCALELES%LKMAX ,   #SCALELES%JMAX -   #SCALELES%IMAX .   #SCALELES%QA /   #SCALELES%KMAX 0   #SCALELES%JA 1   #SCALELES%IA 2   #SCALELES%KA 3   #SCALELES%DAUGHTER_JA 4   #SCALELES%DAUGHTER_IA 5   #SCALELES%DAUGHTER_KA 6   #SCALELES%PARENT_JA 7   #SCALELES%PARENT_IA 8   #SCALELES%PARENT_KA 9   #COMM_WORLD :   #INTERCOMM_PARENT ;   #INTERCOMM_CHILD <   #CNF_FNAME =                                D                                ,                     D                                -                     D                                .                                                      /                     D                                0                                                      1                                                      2                                                      3                                                      4                         p          p            p                                                                    5                         p          p            p                                                                    6                         p          p            p                                                                    7                         p          p            p                                                                    8                         p          p            p                                                                    9                         p          p            p                                    
  @                               :                     
  @                               ;                     
  @                               <                     
  @                              =                              �   0      fn#fn $   �      b   uapp(MOD_LES_DRIVER     �   @   J  SCALE_PRECISION    )  @   J  SCALE_STDIO    i  @   J  SCALE_PROF    �  H   J  DC_LOG    �  M   J  GTOOL_FILE    >  �       LOGINIT+DC_LOG /   �  @      LOGINIT%PRESENT+DC_LOG=PRESENT )     =      LOGINIT%TRIM+DC_LOG=TRIM (   Q  @   a   LOGINIT%FID_CONF+DC_LOG '   �  @   a   LOGINIT%FID_LOG+DC_LOG &   �  @   a   LOGINIT%MASTER+DC_LOG (     H       FILECLOSEALL+GTOOL_FILE "   Y  r       H_MID+SCALE_STDIO #   �  s       H_LONG+SCALE_STDIO %   >  �       IO_SETUP+SCALE_STDIO 5   �  @      IO_SETUP%PRESENT+SCALE_STDIO=PRESENT /   $  =      IO_SETUP%TRIM+SCALE_STDIO=TRIM /   a  P   a   IO_SETUP%MODELNAME+SCALE_STDIO 8   �  @   a   IO_SETUP%CALL_FROM_LAUNCHER+SCALE_STDIO .   �  P   a   IO_SETUP%FNAME_IN+SCALE_STDIO )   A  z       IO_LOG_SETUP+SCALE_STDIO 3   �  =      IO_LOG_SETUP%TRIM+SCALE_STDIO=TRIM 0   �  @   a   IO_LOG_SETUP%MYRANK+SCALE_STDIO 3   8  @   a   IO_LOG_SETUP%IS_MASTER+SCALE_STDIO (   x  @       IO_FID_CONF+SCALE_STDIO '   �  @       IO_FID_LOG+SCALE_STDIO !   �  @       IO_L+SCALE_STDIO &   8	  H       PROF_SETUP+SCALE_PROF )   �	  o       PROF_SETPREFX+SCALE_PROF 3   �	  =      PROF_SETPREFX%TRIM+SCALE_PROF=TRIM 3   ,
  L   a   PROF_SETPREFX%PREFXNAME+SCALE_PROF )   x
  �       PROF_RAPSTART+SCALE_PROF 3     =      PROF_RAPSTART%TRIM+SCALE_PROF=TRIM 9   M  @      PROF_RAPSTART%PRESENT+SCALE_PROF=PRESENT 6   �  L   a   PROF_RAPSTART%RAPNAME_BASE+SCALE_PROF /   �  @   a   PROF_RAPSTART%LEVEL+SCALE_PROF '     �       PROF_RAPEND+SCALE_PROF 1   �  =      PROF_RAPEND%TRIM+SCALE_PROF=TRIM 7   �  @      PROF_RAPEND%PRESENT+SCALE_PROF=PRESENT 4   *  L   a   PROF_RAPEND%RAPNAME_BASE+SCALE_PROF -   v  @   a   PROF_RAPEND%LEVEL+SCALE_PROF *   �  H       PROF_RAPREPORT+SCALE_PROF    �  �      SCALELES ;   �  @     SCALELES%LKMAX+SCALE_LAND_GRID_INDEX=LKMAX 4     @     SCALELES%JMAX+SCALE_GRID_INDEX=JMAX 4   K  @     SCALELES%IMAX+SCALE_GRID_INDEX=IMAX ,   �  @     SCALELES%QA+SCALE_TRACER=QA 4   �  @     SCALELES%KMAX+SCALE_GRID_INDEX=KMAX 0     @     SCALELES%JA+SCALE_GRID_INDEX=JA 0   K  @     SCALELES%IA+SCALE_GRID_INDEX=IA 0   �  @     SCALELES%KA+SCALE_GRID_INDEX=KA A   �  �     SCALELES%DAUGHTER_JA+SCALE_GRID_NEST=DAUGHTER_JA A   _  �     SCALELES%DAUGHTER_IA+SCALE_GRID_NEST=DAUGHTER_IA A   �  �     SCALELES%DAUGHTER_KA+SCALE_GRID_NEST=DAUGHTER_KA =   �  �     SCALELES%PARENT_JA+SCALE_GRID_NEST=PARENT_JA =     �     SCALELES%PARENT_IA+SCALE_GRID_NEST=PARENT_IA =   �  �     SCALELES%PARENT_KA+SCALE_GRID_NEST=PARENT_KA $   C  @   a   SCALELES%COMM_WORLD *   �  @   a   SCALELES%INTERCOMM_PARENT )   �  @   a   SCALELES%INTERCOMM_CHILD #     P   a   SCALELES%CNF_FNAME 