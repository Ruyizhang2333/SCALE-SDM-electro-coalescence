	  �  @   k820309              15.0        ��`                                                                                                           
       common/scale_calendar.F90 SCALE_CALENDAR              CALENDAR_SETUP CALENDAR_GETDAYOFYEAR CALENDAR_DATE2DAYSEC CALENDAR_DAYSEC2DATE CALENDAR_YMD2ABSDAY CALENDAR_HMS2ABSSEC CALENDAR_ADJUST_DAYSEC CALENDAR_COMBINE_DAYSEC CALENDAR_UNIT2SEC CALENDAR_DATE2CHAR I_YEAR I_MONTH I_DAY I_HOUR I_MIN I_SEC                      @                             
                            @                             
                                                                                                                                                                                                                             #         @                                                        #         @                                                       #DAYOFYEAR    #IYEAR 	             D                                     
                 
  @                               	           #         @                                   
                    #ABSDAY    #ABSSEC    #YMDHMS    #SUBSEC    #OFFSET_YEAR              D @                                                     D @                                   
                 
  @                                                      p          p            p                                    
  @                                   
                
  @                                          #         @                                                       #YMDHMS    #SUBSEC    #ABSDAY    #ABSSEC    #OFFSET_YEAR              D @                                                       p          p            p                                    D @                                   
                 
  @                                                    
  @                                   
                
  @                                          #         @                                                     #CALENDAR_YMD2ABSDAY%INT    #CALENDAR_YMD2ABSDAY%MOD    #CALENDAR_YMD2ABSDAY%REAL    #ABSDAY    #GYEAR    #GMONTH    #GDAY    #OYEAR                  @                                INT               @                                MOD               @            @                   REAL           D                                                       
                                                       
                                                       
                                                       
                                             #         @                                                     #CALENDAR_HMS2ABSSEC%REAL     #ABSSEC !   #HOUR "   #MINUTE #   #SECOND $   #SUBSEC %                 @            @                    REAL           D                                !     
                 
  @                               "                     
  @                               #                     
  @                               $                     
                                 %     
      #         @                                   &                   #CALENDAR_ADJUST_DAYSEC%INT '   #CALENDAR_ADJUST_DAYSEC%REAL (   #ABSDAY )   #ABSSEC *                 @                           '     INT               @            @              (     REAL           
D                                 )                      
D                                *     
       %         @                                 +                   
       #CALENDAR_COMBINE_DAYSEC%REAL ,   #ABSDAY -   #ABSSEC .                 @            @              ,     REAL           
  @                               -                     
                                 .     
      #         @                                   /                   #CALENDAR_UNIT2SEC%TRIM 0   #SECOND 1   #VALUE 2   #UNIT 3                 @                           0     TRIM           D                                1     
                 
                                 2     
                
  @                              3                    1 #         @                                   4                    #CHARDATE 5   #YMDHMS 6   #SUBSEC 7   #OFFSET_YEAR 8             D                                5                                      
                                  6                       p          p            p                                    
                                 7     
                
                                  8                                                        9                                                      1                                             :                                                      2                                             ;                                                      3                                             <                                                      4                                             =                                                      5                                             >                                                      6   �   1      fn#fn $   �     b   uapp(SCALE_CALENDAR     �  @   J  SCALE_PRECISION      @   J  SCALE_STDIO #   T  p       DP+SCALE_PRECISION !   �  @       IO_L+SCALE_STDIO '     @       IO_FID_LOG+SCALE_STDIO    D  H       CALENDAR_SETUP &   �  b       CALENDAR_GETDAYOFYEAR 0   �  @   a   CALENDAR_GETDAYOFYEAR%DAYOFYEAR ,   .  @   a   CALENDAR_GETDAYOFYEAR%IYEAR %   n  �       CALENDAR_DATE2DAYSEC ,   �  @   a   CALENDAR_DATE2DAYSEC%ABSDAY ,   7  @   a   CALENDAR_DATE2DAYSEC%ABSSEC ,   w  �   a   CALENDAR_DATE2DAYSEC%YMDHMS ,     @   a   CALENDAR_DATE2DAYSEC%SUBSEC 1   K  @   a   CALENDAR_DATE2DAYSEC%OFFSET_YEAR %   �  �       CALENDAR_DAYSEC2DATE ,     �   a   CALENDAR_DAYSEC2DATE%YMDHMS ,   �  @   a   CALENDAR_DAYSEC2DATE%SUBSEC ,   �  @   a   CALENDAR_DAYSEC2DATE%ABSDAY ,   (  @   a   CALENDAR_DAYSEC2DATE%ABSSEC 1   h  @   a   CALENDAR_DAYSEC2DATE%OFFSET_YEAR $   �  �       CALENDAR_YMD2ABSDAY (   �	  <      CALENDAR_YMD2ABSDAY%INT (   �	  <      CALENDAR_YMD2ABSDAY%MOD )   �	  =      CALENDAR_YMD2ABSDAY%REAL +   5
  @   a   CALENDAR_YMD2ABSDAY%ABSDAY *   u
  @   a   CALENDAR_YMD2ABSDAY%GYEAR +   �
  @   a   CALENDAR_YMD2ABSDAY%GMONTH )   �
  @   a   CALENDAR_YMD2ABSDAY%GDAY *   5  @   a   CALENDAR_YMD2ABSDAY%OYEAR $   u  �       CALENDAR_HMS2ABSSEC )     =      CALENDAR_HMS2ABSSEC%REAL +   R  @   a   CALENDAR_HMS2ABSSEC%ABSSEC )   �  @   a   CALENDAR_HMS2ABSSEC%HOUR +   �  @   a   CALENDAR_HMS2ABSSEC%MINUTE +     @   a   CALENDAR_HMS2ABSSEC%SECOND +   R  @   a   CALENDAR_HMS2ABSSEC%SUBSEC '   �  �       CALENDAR_ADJUST_DAYSEC +   3  <      CALENDAR_ADJUST_DAYSEC%INT ,   o  =      CALENDAR_ADJUST_DAYSEC%REAL .   �  @   a   CALENDAR_ADJUST_DAYSEC%ABSDAY .   �  @   a   CALENDAR_ADJUST_DAYSEC%ABSSEC (   ,  �       CALENDAR_COMBINE_DAYSEC -   �  =      CALENDAR_COMBINE_DAYSEC%REAL /   �  @   a   CALENDAR_COMBINE_DAYSEC%ABSDAY /   3  @   a   CALENDAR_COMBINE_DAYSEC%ABSSEC "   s  �       CALENDAR_UNIT2SEC '   �  =      CALENDAR_UNIT2SEC%TRIM )   5  @   a   CALENDAR_UNIT2SEC%SECOND (   u  @   a   CALENDAR_UNIT2SEC%VALUE '   �  L   a   CALENDAR_UNIT2SEC%UNIT #            CALENDAR_DATE2CHAR ,   �  P   a   CALENDAR_DATE2CHAR%CHARDATE *   �  �   a   CALENDAR_DATE2CHAR%YMDHMS *   d  @   a   CALENDAR_DATE2CHAR%SUBSEC /   �  @   a   CALENDAR_DATE2CHAR%OFFSET_YEAR    �  q       I_YEAR    U  q       I_MONTH    �  q       I_DAY    7  q       I_HOUR    �  q       I_MIN      q       I_SEC 