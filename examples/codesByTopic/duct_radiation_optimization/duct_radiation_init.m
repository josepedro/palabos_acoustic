function duct_radiation_init(radius, mach, reynolds_number, number_processes)
  radius = str2num(radius)
  mach = str2num(mach)
  reynolds_number = str2num(reynolds_number)
  number_processes = str2num(number_processes)

  if reynolds_number <= 5.5148e+03
    command = ['make clean && make && mpirun -n ' num2str(number_processes) ...
    ' duct_radiation_optimization ' num2str(radius) ' ' num2str(mach) ' ' num2str(reynolds_number)];
    disp(command);
    system(command);
  else
    disp('Operation failed')
    disp('Choose Reynolds number < 5.5148e+03');
 end
end
