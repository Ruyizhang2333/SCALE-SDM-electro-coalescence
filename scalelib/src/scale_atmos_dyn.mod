	  ü]     k820309              15.0        ªË`                                                                                                           
       atmos-les/dynamics/scale_atmos_dyn.F90 SCALE_ATMOS_DYN              ATMOS_DYN_SETUP ATMOS_DYN                      @                             
                            @                             
                            @                              
                @          @                              
                            @                             
                @          @                              
                                                                                                                                                                                                                   	                                                       16                                            
                                                                                                                                                                                                                                                                                       1                                                                                                   2                                                                                                   3                                                                                                                                                                                                                                                                              #         @                                                     #PROF_RAPSTART%TRIM    #PROF_RAPSTART%PRESENT    #RAPNAME_BASE    #LEVEL                  @                                TRIM               @                                PRESENT           
                                                     1           
                                                                                                                                                                                                                       #         @                                                     #PROF_RAPEND%TRIM    #PROF_RAPEND%PRESENT    #RAPNAME_BASE     #LEVEL !                 @                                TRIM               @                                PRESENT           
                                                      1           
                                 !                                                        "                                                         #                                                         $                                                         %                                                         &                      D                                  '                      D                                  (            #         @                                   )                  #ATMOS_DYN_SETUP%JA *   #ATMOS_DYN_SETUP%KA +   #ATMOS_DYN_SETUP%IA ,   #ATMOS_DYN_SETUP%SIN -   #DYN_TYPE .   #DENS /   #MOMZ 0   #MOMX 1   #MOMY 2   #RHOT 3   #QTRC 4   #CDZ 5   #CDX 6   #CDY 7   #FDZ 8   #FDX 9   #FDY :   #ENABLE_CORIOLIS ;   #LAT <                                             *                                                      +                                                      ,                          @                           -     SIN           
  @                              .                                    
D @                              /                    
         p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              0                    
         p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              1                    
         p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              2                    
         p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              3                    
         p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              4                    
           p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      5 r        5 r 
     5 r      5 r      5 r                               
  @                              5                    
    p          5 r 
       5 r 
                              
  @                              6                    
    p          5 r        5 r                               
  @                              7                    
    p          5 r        5 r                               
  @                              8                    
     p           5 r 
   n                                       1     5 r 
   n                                      1                                    
  @                              9                    
 !   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              :                    
 "   p           5 r    n                                       1     5 r    n                                      1                                     
                                  ;                    
                                 <                    
 #     p        5 r    p          5 r      5 r        5 r      5 r                      #         @                                   =               C   #ATMOS_DYN%JMAXB >   #ATMOS_DYN%IMAXB ?   #ATMOS_DYN%QA @   #ATMOS_DYN%KMAX A   #ATMOS_DYN%JA B   #ATMOS_DYN%KA C   #ATMOS_DYN%IA D   #ATMOS_DYN%BND_QA E   #ATMOS_DYN%ABS F   #ATMOS_DYN%MAX G   #ATMOS_DYN%REAL H   #DENS I   #MOMZ J   #MOMX K   #MOMY L   #RHOT M   #QTRC N   #DENS_AV O   #MOMZ_AV P   #MOMX_AV Q   #MOMY_AV R   #RHOT_AV S   #QTRC_AV T   #DENS_TP U   #MOMZ_TP V   #MOMX_TP W   #MOMY_TP X   #RHOT_TP Y   #RHOQ_TP Z   #CDZ [   #CDX \   #CDY ]   #FDZ ^   #FDX _   #FDY `   #RCDZ a   #RCDX b   #RCDY c   #RFDZ d   #RFDX e   #RFDY f   #PHI g   #GSQRT h   #J13G i   #J23G j   #J33G k   #MAPF l   #AQ_CV m   #REF_DENS n   #REF_POTT o   #REF_QV p   #REF_PRES q   #ND_COEF r   #ND_COEF_Q s   #ND_ORDER t   #ND_SFC_FACT u   #ND_USE_RS v   #DAMP_DENS w   #DAMP_VELZ x   #DAMP_VELX y   #DAMP_VELY z   #DAMP_POTT {   #DAMP_QTRC |   #DAMP_ALPHA_DENS }   #DAMP_ALPHA_VELZ ~   #DAMP_ALPHA_VELX    #DAMP_ALPHA_VELY    #DAMP_ALPHA_POTT    #DAMP_ALPHA_QTRC    #DIVDMP_COEF    #FLAG_FCT_RHO    #FLAG_FCT_MOMENTUM    #FLAG_FCT_T    #FLAG_FCT_ALONG_STREAM    #USE_AVERAGE    #DTSEC    #DTSEC_ATMOS_DYN    #NSTEP_ATMOS_DYN                                                                                                              >                                                      ?                                                      @                     D                                A                                                      B                                                      C                                                      D                                                        E                          @                           F     ABS               @                           G     MAX               @            @              H     REAL          
D @                              I                    
 $        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              J                    
 %        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              K                    
 &        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              L                    
 '        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              M                    
 (        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D @                              N                    
 )          p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      5 r        5 r 
     5 r      5 r      5 r                               
D                                O                    
 *        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D                                P                    
 +        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D                                Q                    
 ,        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D                                R                    
 -        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D                                S                    
 .        p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
D                                T                    
 /          p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      5 r        5 r 
     5 r      5 r      5 r                               
                                 U                    
 0       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 V                    
 1       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 W                    
 2       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 X                    
 3       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 Y                    
 4       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 Z                    
 5         p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      5 r        5 r 
     5 r      5 r      5 r                               
  @                              [                    
 6   p          5 r 
       5 r 
                              
  @                              \                    
 7   p          5 r        5 r                               
  @                              ]                    
 8   p          5 r        5 r                               
  @                              ^                    
 9   p           5 r 
   n                                       1     5 r 
   n                                      1                                    
  @                              _                    
 :   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              `                    
 ;   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              a                    
 <   p          5 r 
       5 r 
                              
  @                              b                    
 =   p          5 r        5 r                               
  @                              c                    
 >   p          5 r        5 r                               
  @                              d                    
 ?   p           5 r 
   n                                       1     5 r 
   n                                      1                                    
  @                              e                    
 @   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              f                    
 A   p           5 r    n                                       1     5 r    n                                      1                                    
  @                              g                    
 B       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
  @                              h                    
 C         p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      p            5 r 
     5 r      5 r      p                                   
  @                              i                    
 D         p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      p            5 r 
     5 r      5 r      p                                   
  @                              j                    
 E         p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      p            5 r 
     5 r      5 r      p                                    
  @                              k     
               
  @                              l                    
 F         p        p        p        5 r    p        5 r    p          5 r      5 r      p          p            5 r      5 r      p          p                                   
                                 m                    
 G   p          5 r        5 r                               
  @                              n                    
 H       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
  @                              o                    
 I       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
  @                              p                    
 J       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
  @                              q                    
 K       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                                
  @                              r     
                
  @                              s     
                
  @                               t                     
  @                              u     
                
  @                               v                    
                                 w                    
 L       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 x                    
 M       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 y                    
 N       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 z                    
 O       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 {                    
 P       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 |                    
 Q         p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      5 r E       5 r 
     5 r      5 r      5 r E                              
                                 }                    
 R       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                 ~                    
 S       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                                     
 T       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                                     
 U       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                                     
 V       p        5 r    p        5 r 
   p          5 r 
     5 r      5 r        5 r 
     5 r      5 r                               
                                                     
 W         p        5 r    p        5 r    p        5 r 
   p          5 r 
     5 r      5 r      5 r E       5 r 
     5 r      5 r      5 r E                               
  @                                   
                
  @                                                    
  @                                                    
  @                                                    
  @                                                    
                                                       
  @                                   
                
  @                                   
                
  @                                                 ?      fn#fn %   ß   *   b   uapp(SCALE_ATMOS_DYN     	  @   J  SCALE_PRECISION    I  @   J  SCALE_STDIO      @   J  SCALE_PROF !   É  @   j  SCALE_GRID_INDEX    	  @   J  SCALE_INDEX    I  @   j  SCALE_TRACER $     @       IS+SCALE_GRID_INDEX #   É  p       RP+SCALE_PRECISION $   9  r       H_SHORT+SCALE_STDIO $   «  @       KA+SCALE_GRID_INDEX $   ë  @       IA+SCALE_GRID_INDEX $   +  @       JA+SCALE_GRID_INDEX     k  @       QA+SCALE_TRACER &   «  q       ZDIR+SCALE_GRID_INDEX &     q       XDIR+SCALE_GRID_INDEX &     q       YDIR+SCALE_GRID_INDEX !   ş  @       QQA+SCALE_TRACER #   >  p       DP+SCALE_PRECISION !   ®  @       IO_L+SCALE_STDIO '   î  @       IO_FID_LOG+SCALE_STDIO )   .         PROF_RAPSTART+SCALE_PROF 3   Æ  =      PROF_RAPSTART%TRIM+SCALE_PROF=TRIM 9     @      PROF_RAPSTART%PRESENT+SCALE_PROF=PRESENT 6   C  L   a   PROF_RAPSTART%RAPNAME_BASE+SCALE_PROF /     @   a   PROF_RAPSTART%LEVEL+SCALE_PROF !   Ï  @       QQS+SCALE_TRACER !   	  @       QQE+SCALE_TRACER "   O	  @       I_QV+SCALE_TRACER '   	         PROF_RAPEND+SCALE_PROF 1   #
  =      PROF_RAPEND%TRIM+SCALE_PROF=TRIM 7   `
  @      PROF_RAPEND%PRESENT+SCALE_PROF=PRESENT 4    
  L   a   PROF_RAPEND%RAPNAME_BASE+SCALE_PROF -   ì
  @   a   PROF_RAPEND%LEVEL+SCALE_PROF $   ,  @       JS+SCALE_GRID_INDEX $   l  @       JE+SCALE_GRID_INDEX $   ¬  @       IE+SCALE_GRID_INDEX $   ì  @       KS+SCALE_GRID_INDEX $   ,  @       KE+SCALE_GRID_INDEX (   l  @       JBLOCK+SCALE_GRID_INDEX (   ¬  @       IBLOCK+SCALE_GRID_INDEX     ì  G      ATMOS_DYN_SETUP 7   3  @     ATMOS_DYN_SETUP%JA+SCALE_GRID_INDEX=JA 7   s  @     ATMOS_DYN_SETUP%KA+SCALE_GRID_INDEX=KA 7   ³  @     ATMOS_DYN_SETUP%IA+SCALE_GRID_INDEX=IA $   ó  <      ATMOS_DYN_SETUP%SIN )   /  P   a   ATMOS_DYN_SETUP%DYN_TYPE %       a   ATMOS_DYN_SETUP%DENS %       a   ATMOS_DYN_SETUP%MOMZ %   §    a   ATMOS_DYN_SETUP%MOMX %   »    a   ATMOS_DYN_SETUP%MOMY %   Ï    a   ATMOS_DYN_SETUP%RHOT %   ã  T  a   ATMOS_DYN_SETUP%QTRC $   7     a   ATMOS_DYN_SETUP%CDZ $   Ë     a   ATMOS_DYN_SETUP%CDX $   _     a   ATMOS_DYN_SETUP%CDY $   ó    a   ATMOS_DYN_SETUP%FDZ $   ù    a   ATMOS_DYN_SETUP%FDX $   ÿ    a   ATMOS_DYN_SETUP%FDY 0     @   a   ATMOS_DYN_SETUP%ENABLE_CORIOLIS $   E  Ô   a   ATMOS_DYN_SETUP%LAT      ò      ATMOS_DYN 7   !  @     ATMOS_DYN%JMAXB+SCALE_GRID_INDEX=JMAXB 7   K!  @     ATMOS_DYN%IMAXB+SCALE_GRID_INDEX=IMAXB -   !  @     ATMOS_DYN%QA+SCALE_TRACER=QA 5   Ë!  @     ATMOS_DYN%KMAX+SCALE_GRID_INDEX=KMAX 1   "  @     ATMOS_DYN%JA+SCALE_GRID_INDEX=JA 1   K"  @     ATMOS_DYN%KA+SCALE_GRID_INDEX=KA 1   "  @     ATMOS_DYN%IA+SCALE_GRID_INDEX=IA 6   Ë"  @     ATMOS_DYN%BND_QA+SCALE_ATMOS_BOUNDARY    #  <      ATMOS_DYN%ABS    G#  <      ATMOS_DYN%MAX    #  =      ATMOS_DYN%REAL    À#    a   ATMOS_DYN%DENS    Ô$    a   ATMOS_DYN%MOMZ    è%    a   ATMOS_DYN%MOMX    ü&    a   ATMOS_DYN%MOMY    (    a   ATMOS_DYN%RHOT    $)  T  a   ATMOS_DYN%QTRC "   x*    a   ATMOS_DYN%DENS_AV "   +    a   ATMOS_DYN%MOMZ_AV "    ,    a   ATMOS_DYN%MOMX_AV "   ´-    a   ATMOS_DYN%MOMY_AV "   È.    a   ATMOS_DYN%RHOT_AV "   Ü/  T  a   ATMOS_DYN%QTRC_AV "   01    a   ATMOS_DYN%DENS_TP "   D2    a   ATMOS_DYN%MOMZ_TP "   X3    a   ATMOS_DYN%MOMX_TP "   l4    a   ATMOS_DYN%MOMY_TP "   5    a   ATMOS_DYN%RHOT_TP "   6  T  a   ATMOS_DYN%RHOQ_TP    è7     a   ATMOS_DYN%CDZ    |8     a   ATMOS_DYN%CDX    9     a   ATMOS_DYN%CDY    ¤9    a   ATMOS_DYN%FDZ    ª:    a   ATMOS_DYN%FDX    °;    a   ATMOS_DYN%FDY    ¶<     a   ATMOS_DYN%RCDZ    J=     a   ATMOS_DYN%RCDX    Ş=     a   ATMOS_DYN%RCDY    r>    a   ATMOS_DYN%RFDZ    x?    a   ATMOS_DYN%RFDX    ~@    a   ATMOS_DYN%RFDY    A    a   ATMOS_DYN%PHI     B  T  a   ATMOS_DYN%GSQRT    ìC  T  a   ATMOS_DYN%J13G    @E  T  a   ATMOS_DYN%J23G    F  @   a   ATMOS_DYN%J33G    ÔF  T  a   ATMOS_DYN%MAPF     (H     a   ATMOS_DYN%AQ_CV #   ¼H    a   ATMOS_DYN%REF_DENS #   ĞI    a   ATMOS_DYN%REF_POTT !   äJ    a   ATMOS_DYN%REF_QV #   øK    a   ATMOS_DYN%REF_PRES "   M  @   a   ATMOS_DYN%ND_COEF $   LM  @   a   ATMOS_DYN%ND_COEF_Q #   M  @   a   ATMOS_DYN%ND_ORDER &   ÌM  @   a   ATMOS_DYN%ND_SFC_FACT $   N  @   a   ATMOS_DYN%ND_USE_RS $   LN    a   ATMOS_DYN%DAMP_DENS $   `O    a   ATMOS_DYN%DAMP_VELZ $   tP    a   ATMOS_DYN%DAMP_VELX $   Q    a   ATMOS_DYN%DAMP_VELY $   R    a   ATMOS_DYN%DAMP_POTT $   °S  T  a   ATMOS_DYN%DAMP_QTRC *   U    a   ATMOS_DYN%DAMP_ALPHA_DENS *   V    a   ATMOS_DYN%DAMP_ALPHA_VELZ *   ,W    a   ATMOS_DYN%DAMP_ALPHA_VELX *   @X    a   ATMOS_DYN%DAMP_ALPHA_VELY *   TY    a   ATMOS_DYN%DAMP_ALPHA_POTT *   hZ  T  a   ATMOS_DYN%DAMP_ALPHA_QTRC &   ¼[  @   a   ATMOS_DYN%DIVDMP_COEF '   ü[  @   a   ATMOS_DYN%FLAG_FCT_RHO ,   <\  @   a   ATMOS_DYN%FLAG_FCT_MOMENTUM %   |\  @   a   ATMOS_DYN%FLAG_FCT_T 0   ¼\  @   a   ATMOS_DYN%FLAG_FCT_ALONG_STREAM &   ü\  @   a   ATMOS_DYN%USE_AVERAGE     <]  @   a   ATMOS_DYN%DTSEC *   |]  @   a   ATMOS_DYN%DTSEC_ATMOS_DYN *   ¼]  @   a   ATMOS_DYN%NSTEP_ATMOS_DYN 