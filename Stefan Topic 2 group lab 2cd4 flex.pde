TITLE 'IVP numeric Solver' { the problem identification }
COORDINATES cartesian1 { coordinate system, 1D,2D,3D, etc }
VARIABLES { system variables }
theta(threshold=1e-3)
omega(threshold=1e-3)

SELECT { method controls }
ngrid = 1 !low FEM density speeds up solution with no downside when not using FEM
DEFINITIONS { parameter definitions }
alpha =(2560*(cos(2*theta) + 353/8)*((sqrt(16*cos(theta)^2 + 345) + 4*cos(theta))*cos(theta)*sin(theta)*omega^2 - (625*sqrt(16*cos(theta)^2 + 345))/4)*sqrt(353 + 8*cos(2*theta)) + 5120*(sqrt(16*cos(theta)^2 + 345) + 4*cos(theta))*omega^2*(cos(2*theta)^2 + (353*cos(2*theta))/4 + 1)*sin(theta))/(2560*(cos(2*theta) + 353/8)*((cos(theta)^2*sqrt(16*cos(theta)^2 + 345) + 4*cos(theta)^3 - (441*sqrt(16*cos(theta)^2 + 345))/320 - 4*cos(theta))*sqrt(353 + 8*cos(2*theta)) - 2*(sqrt(16*cos(theta)^2 + 345) + 4*cos(theta))*sin(2*theta)*sin(theta)))
INITIAL VALUES
theta = 0
omega = 0

EQUATIONS { PDE's, one for each variable }


theta: dt(theta) = omega
omega: dt(omega) = alpha
! CONSTRAINTS { Integral constraints }

BOUNDARIES { The domain definition }
REGION 1 START(0) LINE TO (1)
TIME 0 TO 10 halt(theta>2*pi) { if time dependent }
!MONITORS { show progress }

PLOTS { save result displays }
for t= 0 by endtime/50 to endtime
history(theta) at(0)
history(omega) at(0)
history(alpha) at(0)
! history(a) at (0)
!history(s) at (0,0) vs v

summary
report t as "time to stop"

END