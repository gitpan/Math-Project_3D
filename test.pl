#!/usr/bin/perl

# Tests for Math::Project_3D and Math::Project_3D::Function
# (c) 2002 Steffen Mueller, all rights reserved

# Forgive me, this is an ugly script.

use strict;
use warnings;

use Test::More tests => 103;

use Math::Project_3D;

ok(1, "Module loaded.");

# Test Math::Project_3D->new_function and MP3D::Function

# Closure syntax first.
{
   # This is a sphere
   my $radius = 5;

   my $x_component = sub {
      my $theta = shift;
      my $phi   = shift;
      return $radius * sin($theta) * cos($phi);
   };
   my $y_component = sub {
      my $theta = shift;
      my $phi   = shift;
      return $radius * sin($theta) * sin($phi);
   };
   my $z_component = sub {
      my $theta = shift;
      return $radius * cos($theta);
   };
 
   my $function = Math::Project_3D->new_function(
     $x_component, $y_component, $z_component         
   );

   ok(ref $function eq 'CODE', "Math::Project_3D->new_function(coderef,...) returned coderef");

   foreach my $u (-2..2) {
      foreach my $v (-2..2) {
         my ($x, $y, $z) = $function->($u, $v);
         ok(
             $x == 5*sin($u)*cos($v) &&
             $y == 5*sin($u)*sin($v) &&
             $z == 5*cos($u),
             "Function returned correct value for ($u, $v)"
         );
      }
   }
}

# Now the plaintext syntax
my $function = Math::Project_3D->new_function(
  'u,v', '$u', 'sin($v)', 'cos($v)' # Looks like a screw
);

ok(ref $function eq 'CODE', "Math::Project_3D->new_function(plaintext) returned coderef");

foreach my $u (-2..2) {
    foreach my $v (-2..2) {
        my ($x, $y, $z) = $function->($u, $v);
        ok(
            $x == $u      &&
            $y == sin($v) &&
            $z == cos($v),
            "Function returned correct value for ($u, $v)"
        );
    }
}

# Now the mixed syntax
$function = Math::Project_3D->new_function(
  'u,v', '$u', sub { sin($_[0]) }, 'cos($v)' # Looks like a screw
);

ok(ref $function eq 'CODE', "Math::Project_3D->new_function(mixed) returned coderef");

foreach my $u (-2..2) {
    foreach my $v (-2..2) {
        my ($x, $y, $z) = $function->($u, $v);

        ok(
            $x == $u      &&
            $y == sin($u) &&
            $z == cos($v),
            "Function returned correct value for ($u, $v)"
        );
    }
}


# Test MP3D->new(), _require_attributes, _make_matrix
# and new_function() on an instance

my $projection = Math::Project_3D->new(
   plane_basis_vector => [ 0, 0, 0 ],
   plane_direction1   => [ .4, 1, 0 ],
   plane_direction2   => [ .4, 0, 1 ],
);

ok(ref $projection eq 'Math::Project_3D', "Math::Project_3D->new returned object");

$function = $projection->new_function('u,v,w', '$u', '$v', '$w');
ok(ref $function eq 'CODE', "MP3D instance->new returned coderef");

# Test get_function and set_function on an instance
ok(ref $projection->get_function() eq 'CODE', "MP3D instance->get_function returned coderef");

# Get a fresh object without a function
$projection = Math::Project_3D->new(
   plane_basis_vector => [ 0, 0, 0 ],
   plane_direction1   => [ .4, 1, 0 ],
   plane_direction2   => [ .4, 0, 1 ],
);
$projection->set_function($function);
ok(ref $projection->get_function() eq 'CODE', "MP3D instance->get_function returned coderef after set_function");


# Now test some projection with this.

# Get a set of data points.

my @result_set;
my $push = 0;
while (<DATA>) {
   chomp;
   last if /^end set\s*$/i and $push;
   push @result_set, split /\s+/ if $push;
   $push = 1 if /^start set\s*$/i;
}

# Above data corresponds to these parameters:
#print join "\n", @result_set;

foreach ( my $u = -1.1; $u <= 1; $u += 0.8 ) {
   foreach ( my $v = 0; $v <= 1; $v += 1 ) {
      foreach ( my $w = 1; $w <= 3; $w += 0.9 ) {
         my ($x, $y, $distance) = $projection->project($u,$v,$w);
         my ($x_c, $y_c, $d_c) = splice @result_set, 0, 3;

         # We're a bit inaccurate. :(
         my $correct = 0;
         $correct = 1 if
            $x        >= $x_c-10E-15 and
            $x        >= $x_c-10E-15 and
            $y        >= $y_c-10E-15 and
            $y        <= $y_c+10E-15 and
            $distance <= $d_c+10E-15 and
            $distance <= $d_c+10E-15;

         ok( $correct,
             "Projected correctly using project() for ($u, $v, $w)."
           );
      }
   }
}


# Now we test project_list()

# Get fresh test data set
$push = 0;
@result_set = ();
while (<DATA>) {
   chomp;
   last if /^end set\s*$/i and $push;
   push @result_set, split /\s+/ if $push;
   $push = 1 if /^start set\s*$/i;
}


# Generate test data
my @array_refs;
foreach ( my $u = -10; $u <= 10; $u += 9 ) {
   foreach ( my $v = 0; $v <= 1; $v += 1 ) {
      foreach ( my $w = 100000; $w <= 1000000; $w += 500000 ) {
         push @array_refs, [$u,$v,$w];
      }
   }
}

my $result_matrix = $projection->project_list(@array_refs);


my $correct = 0;
foreach my $no (1..@array_refs) {
   my ($x_c, $y_c, $d_c) = splice @result_set, 0, 3;

   # == it is not, but eq?????
   $correct++ if
      $result_matrix->element($no, 1) eq $x_c and
      $result_matrix->element($no, 2) eq $y_c and
      $result_matrix->element($no, 3) eq $d_c;
}

ok( $correct == scalar(@array_refs),
    "Projected correctly using project_list()."
  );


# Now we test project_range_callback()

# Get fresh test data set
$push = 0;
@result_set = ();
while (<DATA>) {
   chomp;
   last if /^end set\s*$/i and $push;
   push @result_set, split /\s+/ if $push;
   $push = 1 if /^start set\s*$/i;
}

$projection = Math::Project_3D->new(
  plane_basis_vector => [0,  0, 0],
  plane_direction1 => [.4, 1, 0],
  plane_direction2 => [.4, 0, 1],
  projection_vector  => [1,  1, 1], # defaults to normal of the plane
);

$projection->new_function(
   't,u,v,w,x,y,z', '$t+$x+$w', '$u+$y', '$v+$z',
);

my $sad = 0;
my $okay = 0;
$projection->project_range_callback(
   sub { $sad++; $okay++ unless grep {$_ == shift @result_set} @_ },
   [0, 5, 2.5],
   [-1,1, .69],
   [-1,-1,-1],
   [0,0,-10000],
   [4, 0, -1],
   [1],
   [-1],
);

ok(
    $okay == 45,
    "Projected correctly using project_range_callback ($okay okay of 45)"
  );




__DATA__

The following are the result of the separated project() method
calls.

start set
-0.454545454545455 0.545454545454545 -1.13636363636364
-0.563636363636364 1.33636363636364 -1.40909090909091
-0.672727272727273 2.12727272727273 -1.68181818181818
0.424242424242424 0.424242424242424 -1.43939393939394
0.315151515151515 1.21515151515152 -1.71212121212121
0.206060606060606 2.00606060606061 -1.98484848484848
-0.212121212121212 0.787878787878788 -0.53030303030303
-0.321212121212121 1.57878787878788 -0.803030303030303
-0.43030303030303 2.36969696969697 -1.07575757575758
0.666666666666667 0.666666666666667 -0.833333333333333
0.557575757575758 1.45757575757576 -1.10606060606061
0.448484848484848 2.24848484848485 -1.37878787878788
0.0303030303030303 1.03030303030303 0.0757575757575758
-0.0787878787878788 1.82121212121212 -0.196969696969697
-0.187878787878788 2.61212121212121 -0.46969696969697
0.909090909090909 0.909090909090909 -0.227272727272727
0.8 1.7 -0.5
0.690909090909091 2.49090909090909 -0.772727272727273
end set

This is for the test of the project_list() method
start set
-12124.2424242424 87875.7575757576 -30310.6060606061
-72730.303030303 527269.696969697 -181825.757575758
-12123.3636363636 87875.6363636364 -30310.9090909091
-72729.4242424242 527269.575757576 -181826.060606061
-12121.5151515152 87878.4848484848 -30303.7878787879
-72727.5757575758 527272.424242424 -181818.939393939
-12120.6363636364 87878.3636363636 -30304.0909090909
-72726.696969697 527272.303030303 -181819.242424242
-12118.7878787879 87881.2121212121 -30296.9696969697
-72724.8484848485 527275.151515151 -181812.121212121
-12117.9090909091 87881.0909090909 -30297.2727272727
-72723.9696969697 527275.03030303 -181812.424242424
end set

This is for the test of the project_range_callback() method
start set
-24 -26 24
-36.5 -38.5 36.5
-49 -51 49
-21.93 -24.62 22.62
-34.43 -37.12 35.12
-46.93 -49.62 47.62
-19.86 -23.24 21.24
-32.36 -35.74 33.74
-44.86 -48.24 46.24
-19 -21 19
-31.5 -33.5 31.5
-44 -46 44
-16.93 -19.62 17.62
-29.43 -32.12 30.12
-41.93 -44.62 42.62
-14.86 -18.24 16.24
-27.36 -30.74 28.74
-39.86 -43.24 41.24
-14 -16 14
-26.5 -28.5 26.5
-39 -41 39
-11.93 -14.62 12.62
-24.43 -27.12 25.12
-36.93 -39.62 37.62
-9.86 -13.24 11.24
-22.36 -25.74 23.74
-34.86 -38.24 36.24
-9 -11 9
-21.5 -23.5 21.5
-34 -36 34
-6.93 -9.62 7.62
-19.43 -22.12 20.12
-31.93 -34.62 32.62
-4.86 -8.24 6.24
-17.36 -20.74 18.74
-29.86 -33.24 31.24
-4 -6 4
-16.5 -18.5 16.5
-29 -31 29
-1.93 -4.62 2.62
-14.43 -17.12 15.12
-26.93 -29.62 27.62
0.139999999999999 -3.24 1.24
-12.36 -15.74 13.74
-24.86 -28.24 26.24
end set
