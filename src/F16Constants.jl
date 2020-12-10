# Define Constants for the Aircraft Model

g = 32.17;          # gravity, ft/s^2 
m = 636.94;         # mass, slugs 
B = 30.0;           # span, ft 
S = 300.0;          # planform area, ft^2 
cbar = 11.32;       # mean aero chord, ft 
xcgr = 0.35;        # reference center of gravity as a fraction of cbar 
xcg  = 0.30;        # center of gravity as a fraction of cbar. 
Heng = 0.0;         # turbine momentum along roll axis. 

# NasaData -- translated via eq. 2.4-6 on pg 80 of Stevens and Lewis
Jy  = 55814.0;       # slug-ft^2  
Jxz = 982.0;         # slug-ft^2      
Jz  = 63100.0;       # slug-ft^2 
Jx  = 9496.0;        # slug-ft^2 