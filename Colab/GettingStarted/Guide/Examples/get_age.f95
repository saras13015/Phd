
program get_age
    real :: year, age
    print *, 'What year were you born?'
    read *, year
    age = 2022 - year
    print *, 'Your age is', age
end program get_age
