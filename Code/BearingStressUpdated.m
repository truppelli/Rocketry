% A quantification of the safety factor associated with bearing stress of a
% bolt and the hole that it is screwed into. This is calculated by the
% force applied with respect to area.

clear; 
clc;
f = input('Enter f, where f is the major diameter of the bolt: ');
i = input('Enter i, where i is the inside diameter (minor diameter) of the nut: ');
A=(pi/4).*((f.^2)-(i.^2));
A = abs(A);
disp(A);

P = input('Enter pressure on bolt:');
p = input('Enter pitch of threads: ');
t = input('Enter nut thickness: ');
B=(P./A).*(p./t);
B

T = input('Enter maximum stress before shearing: ')
S=(T./B) % Safety factor calculated

% More information on this analysis can be seen at the link below on pg. 3.
% http://portal.ku.edu.tr/~cbasdogan/Courses/MDesign/course_notes/ScrewStresses.pdf
