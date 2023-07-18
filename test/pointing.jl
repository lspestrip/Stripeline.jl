using Test
using Stripeline

@test [ 0.47450521,-0.4518214 , 0.75544836] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.32727199,-0.50014456, 0.80171595] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.5135005 ,-0.35186432, 0.78262937] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.30988117,-0.42823676, 0.84887392] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.43026602,-0.19561321, 0.88125287] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.37394341,-0.22131343, 0.90065903] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.40920852,-0.45817422, 0.78906576] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.24799695,-0.50279789, 0.82806509] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.47436091,-0.29370891, 0.82988963] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.21422319,-0.3869942 , 0.89685223] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.37522422,-0.14169017, 0.91604076] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.27635352,-0.18526574, 0.94302987] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.36644401,-0.36549968, 0.85564523] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.11400823,-0.43072704, 0.89525211] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.39218549,-0.3009707 , 0.86925668] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.09760689,-0.38734266, 0.91675436] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.23394008,-0.08607117, 0.96843368] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.19176304,-0.10394652, 0.97592113] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.09424787,-0.20455465, 0.97430731] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.00538428,-0.22515833, 0.97430731] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.11583003,-0.14789011, 0.98219749] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [-0.0267358 ,-0.18593891, 0.98219749] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.05702132,-0.08734616, 0.99454473] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.00329291,-0.104259  , 0.99454473] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(45.0)
    )
)) 
            
@test [ 0.48905063,-0.0160092 , 0.87210847] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.33327929,-0.28468507, 0.89882108] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.52173219,0.08710489,0.84865085] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.31130632,-0.31937334, 0.89503578] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.43152168,0.03156518,0.90155015] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.37457005,-0.10740928, 0.92095631] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(1.5),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.42375394,-0.00965638, 0.90572587] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.25400424,-0.28203175, 0.92517022] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.47857304,0.15823242,0.86367259] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.21435008,-0.35346   , 0.91056031] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.37560836,0.09150432,0.92225014] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.276411  ,-0.1501393 , 0.94923925] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(5.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.3746757 ,0.10074025,0.92166671] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.11543339,-0.31688306, 0.94141397] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.39639762,0.16549421,0.90303964] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(20.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.09773378,-0.35311153, 0.93046244] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(20.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.23403818,0.03461074,0.97161115] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.19179246,-0.06768737, 0.9790986 ] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(10.0),
        wheel3ang_0_rad = deg2rad(20.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.09550354,0.04050661,0.99460459] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.00601092,-0.10356438, 0.99460459] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(1.5),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.11621417,0.09770427,0.98840687] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(10.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [-0.02667833,-0.14946613, 0.98840687] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(5.0),
        forkang_rad = deg2rad(10.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [0.05711942,0.03588573,0.9977222 ] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(1.5),
        zVAXang_rad = deg2rad(5.0),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 
            
@test [ 0.00332233,-0.06737489, 0.9977222 ] ≈ rotate_zaxis(telescopetoground(0.0,deg2rad(20.0),deg2rad(20.0),
    TelescopeAngles(
        wheel2ang_0_rad = deg2rad(20.0),
        wheel3ang_0_rad = deg2rad(10.0),
        forkang_rad = deg2rad(5.0),
        zVAXang_rad = deg2rad(1.5),
        ωVAXang_rad = deg2rad(135.0)
    )
)) 