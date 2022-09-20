using Test
using Stripeline

function angletomatrix(wheelanglesfn, time_s, config_ang::configuration_angles_ST)
    rotationmatrix_normalized(camtoground(wheelanglesfn, time_s, config_ang))    
end

@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(rollang_rad=deg2rad(45))),
				[0.2218884684027577 -0.928995249589305 0.29619813272602386; 0.9446039478901318 0.28014092350145736 0.17101007166283433; -0.24184476264797528 0.24184476264797522 0.9396926207859084])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(panang_rad=deg2rad(25))),
				[0.6123724356957945 -0.49999999999999994 0.6123724356957945; 0.3535533905932737 0.8660254037844387 0.3535533905932737; -0.7071067811865476 0.0 0.7071067811865476])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(tiltang_rad=deg2rad(10))),
				[0.8137976813493738 -0.44096961052988237 0.37852230636979245; 0.46984631039295416 0.8825641192593856 0.018028311236297258; -0.3420201433256687 0.16317591116653482 0.9254165783983234])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(rollang_rad=deg2rad(68),tiltang_rad=deg2rad(91),panang_rad=deg2rad(137))),
				[-0.0666774003732695 0.9788974328538561 0.19316816567849124; 0.23777810767845423 -0.17243215916499807 0.9558915848539259; 0.9690282223899963 0.10966752681549406 -0.22126305108141475])
@test isapprox(angletomatrix(_ -> (0, deg2rad(20.0), deg2rad(-30.0)), 0, configuration_angles_ST(wheel2ang_0_rad=deg2rad(48.),wheel3ang_0_rad=deg2rad(-30.),forkang_rad=deg2rad(52.),zVAXang_rad=deg2rad(42.),omegaVAXang_rad=deg2rad(73.),rollang_rad=deg2rad(26),panang_rad=deg2rad(33),tiltang_rad=deg2rad(79))),
				[0.807362184492066 0.3904125522063155 0.4424300420766546; 0.28373613930184904 -0.9143042378965736 0.28903557535093766; 0.517358779094625 -0.10782300146126961 -0.8489488170965328])
