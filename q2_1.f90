program q2
  implicit none
  real(8) :: eta(0:1000), u(0:1000)
  real(8) :: eta_1(0:1000), eta_2(0:1000), eta_3(0:1000)
  real(8) :: u_1(0:1000), u_2(0:1000), u_3(0:1000)
  real(8) :: dt = 0.01_8, dx = 0.1_8, h = 1.0_8, g = 9.8_8
  real(8), parameter :: pi = Atan(1.0_8) * 4.0_8
  integer :: n = 10000, m = 1000, i, j ! n,i: time; m,j: space

  ! initial condition
  do j = 0, m
     eta(j) = exp(-100.0_8 * (real(j,8) / real(m,8) - 0.5_8)**2)
     u(j) = 0.0_8
  end do

  ! Output
  do j = 0, m
     write(*,*) dx * j, eta(j), u(j)
  end do
  write(*,*)

  ! Time evolution using Runge-Kutta method
  do i = 1, n
     do j = 1, m-1
        eta_1(j) = eta(j) - dt * h * (u(j+1) - u(j-1)) / (4.0_8 * dx)
        u_1(j) = u(j) - dt * g * (eta(j+1) - eta(j-1)) / (4.0_8 * dx)
     end do
     eta_1(0) = eta(1)
     eta_1(m) = eta(m-1)
     u_1(0) = 0.0_8
     u_1(m) = 0.0_8
     do j = 1, m-1
        eta_2(j) = eta(j) - dt * h * (u_1(j+1) - u_1(j-1)) / (4.0_8 * dx)
        u_2(j) = u(j) - dt * g * (eta_1(j+1) - eta_1(j-1)) / (4.0_8 * dx)
     end do
     eta_2(0) = eta_1(1)
     eta_2(m) = eta_1(m-1)
     u_2(0) = 0.0_8
     u_2(m) = 0.0_8
     do j = 1, m-1
        eta_3(j) = eta(j) - dt * h * (u_2(j+1) - u_2(j-1)) / (2.0_8 * dx)
        u_3(j) = u(j) - dt * g * (eta_2(j+1) - eta_2(j-1)) / (2.0_8 * dx)
     end do
     eta_3(0) = eta_2(1)
     eta_3(m) = eta_2(m-1)
     u_3(0) = 0.0_8
     u_3(m) = 0.0_8
     do j = 1, m-1
        eta(j) = eta(j) - dt * h * ((u(j+1) - u(j-1)) + (u_3(j+1) - u_3(j-1)) &
            & + 2.0_8 * (u_1(j+1) - u_1(j-1) + u_2(j+1) - u_2(j-1))) / (12.0_8 * dx)
        u(j) = u(j) - dt * g * ((eta(j+1) - eta(j-1)) + (eta_3(j+1) - eta_3(j-1)) &
            & + 2.0_8 * (eta_1(j+1) - eta_1(j-1) + eta_2(j+1) - eta_2(j-1))) / (12.0_8 * dx)
     end do
     eta(0) = eta(1)
     eta(m) = eta(m-1)
     u(0) = 0.0_8
     u(m) = 0.0_8

     ! Output
     if (mod(i,100) == 0) then
        do j = 0, m
           write(*,*) dx * j, eta(j), u(j)
        end do
        write(*,*)
     end if
  end do

  stop
end program q2
