
# See the POD documentation at the end of this
# document for detailed copyright information.
# (c) 2002 Steffen Mueller, all rights reserved.

package Math::Project_3D;

use strict;
use warnings;

use 5.006;

use vars qw/$VERSION/;

$VERSION = 1.004;

use Carp;

use Math::MatrixReal;
use Math::Project_3D::Function;


# class and object method new_function
# 
# Wrapper around Math::Project_3D->new()
# Returns compiled function (anon sub)
# As a class method, it does not have side effects.
# As an object method, it assigns the returned function to
# the object's function attribute.

sub new_function {
   @_ > 1 or croak "No arguments supplied to 'new_function' method.";

   my $self = shift;

   my @components = @_;

   my $function = Math::Project_3D::Function->new(@components);

   if (ref $self eq __PACKAGE__) {
      $self->{function} = $function;
   }

   return $function;
}


# Class and object method/constructor new
# 
# Arguments are used in the object's anon hash.
# Creates new object.
# Returns MP3D instance.

sub new {
   my $proto = shift;
   my $class = ref $proto || $proto;

   # !!!FIXME!!! Insert argument checking here
   my $self  = {
     function => undef,
     @_
   };

   bless $self => $class;

   # Some attributes are required.
   my $missing_attribute = $self->_require_attributes(qw(
                                   plane_basis_vector
                                   plane_direction1
                                   plane_direction2
                           ));
   croak "Required attribute '$missing_attribute' missing."
     if defined $missing_attribute;

   # Transform the vector/matrix specifications ((nested) array refs)
   # into Math::MatrixReal's.
   foreach ( qw(
        plane_basis_vector
        plane_direction1
        plane_direction2
      ) ) {

      # Argument checking, we either want an array ref or a MMR object.
      croak "Invalid argument '$_' passed to new. Should be an array reference."
        if not exists $self->{$_};
      
      next if ref $self->{$_} eq 'Math::MatrixReal';
      croak "Invalid argument '$_' passed to new. Should be an array reference."
        if ref $self->{$_} ne 'ARRAY';

      # Transform into an MMR object.
      $self->{$_} = $self->_make_matrix($self->{$_});
   }

   # Projection vector defaults to the normal vector of the plane,
   # but may also be specified explicitly as array ref or MMR object.
   if ( defined $self->{projection_vector} and
        ref $self->{projection_vector} eq 'ARRAY' ) {

      $self->{projection_vector} = $self->_make_matrix( $self->{projection_vector} );

   } elsif ( not defined $self->{projection_vector} or
             ref $self->{projection_vector} ne 'Math::MatrixReal' ) {

      # Defaults to the normal of the plane.
      $self->{projection_vector} = $self->{plane_direction1}->vector_product(
                                     $self->{plane_direction2}
                                   );
   }

   # Now generate the linear equation system that has to be solved by
   # MMR in order to do project a point onto the plane.
   # 
   # MMR does all the dirty work. (Thanks!) So we just create a matrix.
   # (d, e be the directional vectors, n the vector we project along,
   #  p the basis point. Let i be the solution)
   # 
   # x(t) + n1*i1 = d1*i2 + e1*i3 + p1
   # y(t) + n2*i1 = d2*i2 + e2*i3 + p2
   # z(t) + n3*i1 = d3*i2 + e3*i3 + p3
   # 
   # This is the intersection of the plane and the corresponding orthogonal
   # line through s(t). Another form:
   # 
   # n1*i1 + d1*i2 + e1*i3 = p1 + x(t)
   # n2*i1 + d2*i2 + e2*i3 = p2 + y(t)
   # n3*i1 + d3*i2 + e3*i3 = p3 + z(t)
   # 
   # linear_eq_system will be the matrix of n,d,e.
   # result_vector will be p+s(t) (later)

   my $linear_eq_system = Math::MatrixReal->new_from_cols(
                         [
                           $self->{projection_vector},
                           $self->{plane_direction1},
                           $self->{plane_direction2},
                         ]
   );

   # Now we do the first step towards solving the equation system.
   # This does not have to be repeated for every projection.
   $self->{lr_matrix} = $linear_eq_system->decompose_LR();

   return $self;
}


# Method project
# 
# Does the projection calculations for a single set of function
# parameters.
# Takes the function parameters as argument.
# Returns undef on failure (plane and P + lambda*projection_vector do
# not intersect!). On success:
# Returns the coefficients of the point on the plane, 
# the distance of the source point from the plane in lengths
# of the projection vector.

sub project {
   my $self  = shift;

   # Apply function
   my $point = Math::MatrixReal->new_from_cols(
                 [
                   [ $self->{function}->(@_) ],
                 ]
   );

   # Generate result_vector
   my $result_vector = Math::MatrixReal->new_from_cols(
                         [
                           $self->{plane_basis_vector} + $point
                         ]
   );

   # Solve the l_e_s.
   my ($dimension, $p_vector, undef) = $self->{lr_matrix}->solve_LR($result_vector);
   
   # Did we find a solution?
   return undef if not defined $dimension;

   # $dimension == 0 => one solution
   # $dimension == 1 => straight line
   # $dimension == 2 => plane (impossible :) )
   # ...
   # $dimension == 1 is possible (point and projection_vector part
   # of the plane). Hence: !!!FIXME!!!

   return $p_vector->element(2,1), # coefficient 1
          $p_vector->element(3,1), # coefficient 2
          $p_vector->element(1,1); # distance in lengths of the projection vector
}


# Method project_range_callback
# 
# Does the projection calculations for ranges of function parameters.
# Takes an code reference and a number array refs of ranges as argument.
# Ranges are specified using three values:
# [lower_bound, upper_bound, increment].
# Alternatively, one may specify only one element in any of the array
# references. This will then be a static parameter.
# The number of ranges (array references) corresponds to the number of
# parameters to the vectorial function.
# The callback is called for every set of result values with an array
# reference of parameters and the list of the three result values.

sub project_range_callback {
   my  $self     = shift;

   # Argument checking
   my  $callback = shift;
   ref $callback eq 'CODE' or
     croak "Invalid code reference passed to project_range_callback.";

   my @ranges    = @_;

   croak "Invalid range parameters passed to project_range_callback"
     if grep { ref $_ ne 'ARRAY' } @ranges;


   my $function = $self->get_function();

   # Replace all ranges that consist of a single number with ranges
   # that iterate from that number to itself.
   foreach (@ranges) {
      $_ = [$_->[0], $_->[0], 1] if @$_ == 1;
   }

   # Calculate every range's length and store it in @lengths.
   my @lengths   = map {
                         #      upper  - lower    / increment
                         int( ($_->[1] - $_->[0]) / $_->[2] ),
                       } @ranges;

   # Prepare counters for every range.
   my @counters  = (0) x scalar(@ranges);

   # Calculate the number if iterations needed.
   # It is $_+1 and not $_ because the lengths are the lengths we
   # need for the comparisons inside the long for loop. We save one
   # op in there that way.
   my $iterations = 1;
   $iterations *= ( $_ + 1 ) for @lengths;
   
   # For all possible combinations of parameters...
   for (my $i = 1; $i <= $iterations; $i++) {
      
      # Get current function parameters
      my @params;

      # Get one parameter for every range
      for (my $range_no = 0; $range_no < @ranges; $range_no++) {
         # lower + increment * current_count
         push @params, $ranges[$range_no][0] +
                       $ranges[$range_no][2] * $counters[$range_no];
      }
      
      # Increment outermost range by one. If it got out of bounds,
      # make it 0 and increment the next range, etc.
      my $j = 0;
      while (defined $ranges[$j] and ++$counters[$j] > $lengths[$j]) {
         $counters[$j] = 0;
         $j++;
      }

      # Apply function
      my $point = Math::MatrixReal->new_from_cols(
                    [
                      [ $function->(@params) ],
                    ]
      );

      # Generate result_vector
      my $result_vector = Math::MatrixReal->new_from_cols(
                            [
                              $self->{plane_basis_vector} + $point
                            ]
      );

      # Solve the l_e_s.
      my ($dimension, $p_vector, undef) = $self->{lr_matrix}->solve_LR($result_vector);
   
      # Did we find a solution?
      croak "Could not project $result_vector."
        if not defined $dimension;

      # $dimension == 0 => one solution
      # $dimension == 1 => straight line
      # $dimension == 2 => plane (impossible :) )
      # ...
      # $dimension == 1 is possible (point and projection_vector part
      # of the plane). Hence: !!!FIXME!!!
      
      $callback->(
             $p_vector->element(2,1), # coefficient 1
             $p_vector->element(3,1), # coefficient 2
             $p_vector->element(1,1), # distance in lengths of the projection vector
      );
   }
   
   return();
}


# Method project_list
# 
# Wrapper around project(), therefore slow.
# Takes a list of array refs as argument. The array refs
# are to contain sets of function parameters.
# Calculates every set of parameters' projection and stores the
# three associated values (coefficients and distance coefficient)
# in an n*3 matrix (MMR obj) wheren is the number of points
# projected.
# Currently croaks if any of the points cannot be projected onto
# the plane using the given projection vector. (In R3 -> R2, try
# using the plane's normal vector which is guaranteed not to be
# parallel to the plane.)
# Returns that MMR object.

sub project_list {
   my $self = shift;
   croak "No arguments passed to project_list()."
     if @_ == 0;

   # Create result matrix to hold individual results.
   my $result_matrix = Math::MatrixReal->new(scalar(@_), 3);

   # Keep track of the matrix row
   my $result_no = 1;

   foreach my $array_ref (@_) {
      my ($coeff1, $coeff2, $dist) = $self->project(@$array_ref);

      croak "Vector $result_no cannot be projected."
        if not defined $coeff1;

      # Assign results
      $result_matrix->assign($result_no, 1, $coeff1);
      $result_matrix->assign($result_no, 2, $coeff2);
      $result_matrix->assign($result_no, 3, $dist);
      $result_no++;
   }

   return $result_matrix;
}


# Accessor get_function
# 
# No parameters.
# Returns the current function code ref.

sub get_function {
   my $self = shift;
   return $self->{function};
}


# Accessor set_function
# 
# Takes a code ref as argument.
# Sets the object's function to that code ref.
# Returns the code ref.

sub set_function {
   my $self     = shift;
   my $function = shift;

   ref $function eq 'CODE' or croak "Argument to set_function must be code reference.";

   $self->{function} = $function;

   return $self->{function};
}


# Private method _require_attributes
# 
# Arguments must be a list of attribute names (strings).
# Tests for the existance of those attributes.
# Returns the missing attribute on failure, undef on success.

sub _require_attributes {
   my $self = shift;
   exists $self->{$_} or return $_ foreach @_;
   return undef;
}


# Private method _make_matrix
# 
# Takes a list of array refs as arguments.
# Creates a Math::MatrixReal object from the arrays.
# Returns the Math::MatrixReal object.

sub _make_matrix {
   my $self = shift;
   @_ or croak "No arguments passed to _make_matrix.";

   my $matrix = Math::MatrixReal->new_from_cols( [ @_ ] );

   return $matrix;
}


1;

__END__

=pod

=head1 NAME

Math::Project_3D - Project functions of multiple parameters
from R^3 onto an arbitrary plane

=head1 VERSION

Current version is 1.002. Beta. Use at your own risk.

=head1 SYNOPSIS

  use Math::Project_3D;
  
  my $projection = Math::Project_3D->new( {
    plane_basis_vector => [0,  0, 0],
    plane_direction1   => [.4, 1, 0],
    plane_direction2   => [.4, 0, 1],
    projection_vector  => [1,  1, 1], # defaults to normal of the plane
  } );

  $projection->new_function(
    'u,v', 'sin($u)', 'cos($v)', '$u' 
  );

  ($plane_coeff1, $plane_coeff2, $distance_coeff) =
     $projection->project( $u, $v );

=head1 NOTE

This module is beta software. That means it should work as expected
(by the author) and as documented, but it has not been thoroughly
enough tested in order to qualify as "stable" or "mature" code. The
documentation is still sparse (which is why it works as documented),
but that should be fixed real-soon-now.

I'll start explaining what this module does with some background. Feel
free to skip to L<DESCRIPTION> if you don't feel like vector geometry.

=head1 BACKGROUND

Given a function of three components and of an arbitrary number of
parameters, plus a few vectors, this module creates a projection of
individual points on this vectorial function onto an arbitrary plane
in three dimensions.

The module does this by creating lines from the result of the vectorial
function s(a) = x,y,z and a specified projection vector (which defaults
to the normal vector of the projection plane. The normal vector is defined
as being orthogonal to both directional vectors of the plane or as the
vector product of the two.). Then, using the linear equation solver of
Math::MatrixReal, it calculates the intersection of the line and the plane.

This point of intersection can be expressed as

  basis point of the plane + i2 * d + i3 * e

where i2/i3 are the coefficients that are the solution we got from
solving the linear equation system and d1/d2 are the directional
vectors of the plane. Basically, the equation system looks like this:

   n1*i1 + d1*i2 + e1*i3 = p1 + x(t)
   n2*i1 + d2*i2 + e2*i3 = p2 + y(t)
   n3*i1 + d3*i2 + e3*i3 = p3 + z(t)

where n1/2/3 are the normal vector components. p1/2/3 the basis point
components, t is a vector of function parameters. i the solution.

Now, on the plane, you can express the projected point in terms of
the directional vectors and the calculated coefficients.

=head1 DESCRIPTION

=head2 Methods

=over 4

=item new

C<new> may be used as a class or object method. In the current
implementation both have the same effect.

C<new> returns a new Math::Project_3D object. You need to pass a
number of arguments as a list of key/value pairs:

  plane_basis_vector => [0,  0, 0],
  plane_directional1 => [.4, 1, 0],
  plane_directional2 => [.4, 0, 1],
  projection_vector  => [1,  1, 1], # defaults to normal of the plane
  
plane_basis vector denotes the position of the basis point of the plane
you want to project onto. Vectors are generally passed as array references
that contain the components. Another way to pass vectors would be to pass
Math::MatrixReal objects instead of the array references.
plane_directional1 and plane_directional2 are the vectors that span the
plane. Hence, the projection plane has the form:

  s = plane_basis_vector + coeff1 * plane_dir1 + coeff2 * plane_dir2

The last vector you need to specify at creation time is the vector along
which you want to project the function points onto the plane. You may,
however, omit its specification because it defaults to teh cross-product
of the plane's directional vectors. Hence, all points are orthogonally
projected onto the plane.

=item new_function

This method may be used as a class or object method. It does not have
side effects as a clas method, but as an object method, its results
are applied to the projection object.

new_function returns an anonymous subroutine compiled from component
functions which you need to specify.

For a quick synopsis, you may look at L<Math::Project_3D::Function>.

You may pass a list of component functions that are included in the
compiled vectorial functions in the order they were passed. There are
two possible ways to specify a component function. Either you pass
an subroutine reference which is called with the list of parameters,
or you pass a string containing a valid Perl expression which is then
evaluated. You may mix the two syntaxes at will.

If any one of the component functions is specified as a string,
the first argument to new_function I<must> be a string of parameter
names separated by commas. These parameter names will then be made
availlable to the string component functions as the respective
scalar variables. (eg. 't,u' will mean that the parameters availlable
to the string expressions are $t and $u.)

Due to some black magic in Math::Project_3D::Function, the string
expression syntax is actually faster at run-time because you save
subroutine calls which are a major bottleneck in Perl.
(Like 50 ops or so?) Arguably, the closure syntax is more powerful
because being a closure, it has access to variables I<outside> the
scope of the resulting vectorial function. For a simple-minded
example, you may have a look at the synopsis in
L<Math::Project_3D::Function>. Picture a dynamic radius, etc.

=item get_function

get_function returns the object's current vectorial function.

=item set_function

set_function is the counterpart of get_function.

=item project

The project method can be used to do the projection calculations for
one point of the vectorial function. It expects a list of function
parameters as argument.

On failure (eg. the projection vector is parallel to the plane.),
the method returns undef.

On success, the method returns the coefficients of the projected
point on the plane as well as the distance of the source point from
the plane in lengths of the projection vector. (May be negative for
points "behind" the plane.)

=item project_list

project_list is a wrapper around project and therefore rather slow
because it is doing a lot of extra work.

project_list takes a list of array references as arguments. The
referenced arrays are to contain sets of function parameters.
The method the calculates every set of parameters' projection and
stores the three results (coefficients on the plane and distance
coefficient) in an n*3 matrix (as an MMR object, n is the number of
points projected). The matrix is returned.

Currently, the method croaks if any of the points cannot be
projected onto the plane using the given projection vector.
The normal vector of the plane used as the projection vector should
guarantee valid results for any points.

=item project_range_callback

For projection of a large number of points, this method will probably
be the best bet. Its first argument has to be a callback function that
will be called with the calculated coefficients for every projected point.

All arguments thereafter have to be array references. Every one of
these referenced arrays represents the range of one parameter.
These arrays may either contain one number, which is then treated as a
static parameter, or three numbers of the form:

  [ lower boundary, upper boundary, increment ]

For example, [-1, 1, 0.8] would yield the parameter values
-1, -0.2, 0.6. You may also reverse upper and lower boundary,
given you increment by a negative value: [1, -1, -0.8] yields the
same parameter values but in reversed order. Example:

  $projection->project_range_callback(
    sub {...}, # Do some work with the results
    [ 1, 2,  .5],
    [ 2, 1,  .5],
    [ 0, 10, 4 ],
  );

Will result in the coefficients for the following parameter
sets being calculated:

1,2,0 1.5,2,0 2,2,0 1,1.5,0 1.5,1.5,0 2,1.5,0 1,1,0 1.5,1,0 2,1,0
1,2,4 etc.

croaks if a point cannot be projected. (projection vector parallel
to the plane but not I<in> the plane.)

=back

=head1 AUTHOR

Steffen Mueller, mail at steffen-mueller dot net

=head1 COPYRIGHT

Copyright (c) 2002 Steffen Mueller. All rights reserved.
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 SEE ALSO

L<Math::MatrixReal>

=cut
