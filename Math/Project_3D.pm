
# See the POD documentation at the end of this
# document for detailed copyright information.
# (c) 2002 Steffen Mueller, all rights reserved.

package Math::Project_3D;

use strict;
use warnings;

use vars qw/$VERSION/;

$VERSION = 1.002;

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
    parameters         => 't,u',
    plane_basis_vector => [0,  0, 0],
    plane_directional1 => [.4, 1, 0],
    plane_directional2 => [.4, 0, 1],
    projection_vector  => [1,  1, 1], # defaults to normal of the plane
  } );

  $projection->new_function(
    'u,v', 'sin($u)', 'cos($v)', '$u' 
  );

  ($plane_coeff1, $plane_coeff2, $distance_coeff) =
     $projection->project( $u, $v );

=head1 DESCRIPTION

This module is beta software. That means it should work as expected
(by the author) and as documented, but it has not been thoroughly
enough tested in order to qualify as "stable" or "mature" code. The
documentation is still sparse (which is why it works as documented),
but that should be fixed real-soon-now.

=head1 AUTHOR

Steffen Mueller, mail at steffen-mueller dot net

=head1 COPYRIGHT

Copyright (c) 2002 Steffen Mueller. All rights reserved.
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 SEE ALSO

L<Math::MatrixReal>

=cut
