# Writing an INPUT file

There is no obligation to put any order in the keywords

## MANDATORY KEYWORDS

* BATH bath_radius

* PSEUDO first_shell_radius

* OUTPUT output_file_name

* PATTERN nAtom1 labelAtom1 nAtom2 labelAtom2 ...  

Not that, if you have more than once the Atom1, you will have to place it again at the end (for example, 5 Ir 8 Sr 6 O 5 Ir)

* NPATTERN nPattern

* LATTICE

  a cell_parameter_a
  b cell_parameter_b
  c cell_parameter_c
  alpha cell_angle_alpha
  beta cell_angle_beta
  gamma cell_angle_gamma

* ATOMS label charge neighbours radius_of_the_coordination_sphere ...

* SYMOP
   x,y,z
   ...
   POMYS
   
* GENERATOR
  atom1 x_frac y_frac z_frac
  ...
  ROTARENEG
  
* END_OF_INPUT


## OPTIONNAL KEYWORDS

* TRANSLATE x y z

Crystal MET has some trouble finding specific pattern, sometimes forcing a translation can help finding it

* CENTER  atom1 atom2 ...

Will chose the atom (or the center of the atoms) as the center of the sphere

* X_AXIS atom1 atom2 ...

Will chose the atom (or the center of the atoms) as the x axis

* Y_AXIS  atom1  atom2 ...

Will chose the atom (or the center of the atoms) as the y axis

* Z_AXIS atom1  atom2 ...

Will chose the atom (or the center of the atoms) as the z axis

* NOTINPP atoms not in the first shell

* COLOR  

Shows the bath sorted with colors

* NOCOLOR 

Shows the bath sorted by atoms

* SEEFRAG  

Shows the fragment once found (recommended)

* OPTIMIZATION  

Will use an optimization method to set the total charge to 0 by slightly changing the charges of the outer shell atoms

* EVJEN 
Simili Evjen method to calculate the charge of the atom on the outer shell, note that OPTIMIZATION and EVJEN keywords can not be used together

