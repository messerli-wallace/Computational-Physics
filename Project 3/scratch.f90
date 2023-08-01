DO ball_spin_index = 0, 3000 !Test ball spin 0-300 rad/s with interval of 0.1
  if (test .NE. 0) THEN
    omega = ball_spin_index * (0.1d0)
  else
    omega = 0.0d0
    DO theta_index = 0, 3000 !Test from 0-30 degrees
      current_angle = (theta_index)*0.01d0 ! current angle
  
      !call calculate for current angle
      Call calculate(x,y,z,vx,vy,vz,current_angle,n,omega)
      write(1,*) current_angle,x(n) !write the angle and final distance to track in file
      !check if current last distance is greater than our desired amount
      IF (x(n) .gt. final_x) THEN
        print*, current_angle, x(n) !print the straddling distances and their angle
        print*, current_angle-0.1d0, previous_dist
        exit
      END IF
      previous_dist = x(n) !store distance in case the next angle is the last one
    END DO
  END IF

  if (test .EQ. 0) exit !exit the omega loop if we aren't testing for omega
END DO

Call output(x,y,z,vx,vy,vz,n,current_angle,omega)