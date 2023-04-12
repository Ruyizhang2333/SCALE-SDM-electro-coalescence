#!/bin/bash

cat << EOF > conf/param.region.conf

#################################################
#
# model configuration: process
#
#################################################

&PARAM_PRC
 PRC_NUM_X      = 4,
 PRC_NUM_Y      = 4,
 PRC_PERIODIC_X = .false.,
 PRC_PERIODIC_Y = .false.,
/

#################################################
#
# model configuration: region
#
#################################################

&PARAM_INDEX
 KMAX = 80,
 IMAX = 320,
 JMAX = 320,
/

&PARAM_LAND_INDEX
 LKMAX = 5,
/

&PARAM_URBAN_INDEX
 UKMAX = 4,
/

&PARAM_LAND_GRID
 LDZ = 0.05D0, 0.15D0, 0.30D0, 0.50D0, 1.00D0,
/

&PARAM_URBAN_GRID
 UDZ = 0.05D0, 0.15D0, 0.30D0, 0.50D0,
/

&PARAM_GRID
 DX = 100.D0,
 DY = 100.D0,
 FZ(:) =    100.0000D0,   200.0000D0,   300.0000D0,   400.0000D0,   500.0000D0,
            600.0000D0,   700.0000D0,   800.0000D0,   900.0000D0,  1000.0000D0,
           1100.0000D0,  1200.0000D0,  1300.0000D0,  1400.0000D0,  1500.0000D0,
           1600.0000D0,  1700.0000D0,  1800.0000D0,  1900.0000D0,  2000.0000D0,
           2112.0000D0,  2230.2720D0,  2355.1672D0,  2487.0566D0,  2626.3318D0,
           2773.4062D0,  2928.7170D0,  3092.7251D0,  3265.9177D0,  3448.8091D0,
           3641.9424D0,  3845.8911D0,  4061.2610D0,  4288.6914D0,  4528.8579D0,
           4782.4741D0,  5050.2925D0,  5333.1089D0,  5631.7632D0,  5947.1421D0,
           6280.1821D0,  6631.8726D0,  7003.2573D0,  7395.4399D0,  7809.5845D0,
           8246.9209D0,  8708.7480D0,  9196.4375D0,  9711.4385D0, 10255.2793D0,
          10805.2793D0, 11355.2793D0, 11905.2793D0, 12455.2793D0, 13005.2793D0,
          13555.2793D0, 14105.2793D0, 14655.2793D0, 15205.2793D0, 15755.2793D0,
          16305.2793D0, 16855.2793D0, 17405.2793D0, 17955.2793D0, 18505.2793D0,
          19055.2793D0, 19605.2793D0, 20155.2793D0, 20705.2793D0, 21255.2793D0,
          21784.0234D0, 22290.9844D0, 22775.6543D0, 23237.5488D0, 23676.2051D0,
          24091.1855D0, 24482.0742D0, 24848.4805D0, 25190.0391D0, 25506.4082D0,
 BUFFER_DZ =   0.D0,
 BUFFER_DX = 800.D+3,
 BUFFER_DY = 800.D+3,
 BUFFFACT  =   1.0255605111453550D0,
/

&PARAM_MAPPROJ
 MPRJ_basepoint_lon = 135.220404D0,
 MPRJ_basepoint_lat =  34.653396D0,
/
EOF