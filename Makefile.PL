use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Math::Project_3D',
    'VERSION_FROM'	=> 'Math/Project_3D.pm', # finds $VERSION
    'PREREQ_PM'		=> {Math::MatrixReal => 1.9},
    ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM => 'Math/Project_3D.pm', # retrieve abstract from module
       AUTHOR     => 'Steffen Mueller <project_3d-module@steffen-mueller.net>') : ()),
);
