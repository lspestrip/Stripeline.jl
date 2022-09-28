using Test
using Stripeline

function angletomatrix(wheelanglesfn, time_s, config_ang::configuration_angles_ST)
    rotationmatrix_normalized(camtoground(wheelanglesfn, time_s, config_ang))    
end

@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(rollang_rad=deg2rad(45))),
				[0.2218884684027577 -0.928995249589305 0.29619813272602386; 0.9446039478901318 0.28014092350145736 0.17101007166283433; -0.24184476264797528 0.24184476264797522 0.9396926207859084])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(panang_rad=deg2rad(25))),
				[0.8137976813493738 -0.32797515353481177 0.4797558050656603; 0.46984631039295416 0.8571575464476954 -0.21101039116094453; -0.3420201433256687 0.39713126196710286 0.8516507396391465])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(tiltang_rad=deg2rad(10))),
				[0.7500000000000001 -0.49999999999999994 0.4330127018922193; 0.4330127018922193 0.8660254037844387 0.24999999999999994; -0.49999999999999994 0.0 0.8660254037844387])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(rollang_rad=deg2rad(68),tiltang_rad=deg2rad(91),panang_rad=deg2rad(137))),
				[0.47444246561235864 0.34112613016341686 0.8115031177656669; -0.21412297968452368 -0.8494536195128901 0.4822653811621473; 0.8538475838196693 -0.4025686421173186 -0.3299739262261965])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(wheel2ang_0_rad=deg2rad(48.),wheel3ang_0_rad=deg2rad(-30.),forkang_rad=deg2rad(52.),zVAXang_rad=deg2rad(42.),omegaVAXang_rad=deg2rad(73.),rollang_rad=deg2rad(26),panang_rad=deg2rad(33),tiltang_rad=deg2rad(79))),
				[-0.17570296000187266 0.5751787563280825 0.7989354592928397; 0.45810521448987807 -0.670565309868962 0.5835081641738409; 0.8713599040027968 0.4685206115735181 -0.14567207771914797])
