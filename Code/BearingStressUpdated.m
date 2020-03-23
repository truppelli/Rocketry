% A quantification of the safety factor associated with bearing stress of a
% bolt and the hole that it is screwed into. This is calculated by the
% force applied with respect to area.

clear; 
clc;
f = input('where f is the major diameter of the bolt');
i = input('where i is the inside diameter (minor diameter) of the nut');
A=(pi/4).*((f.^2)-(i.^2));
A = abs(A);
disp(A);

P = input('is pressure on bolt');
p = input('is pitch of threads');
t = input('is nut thickness');
B=(P./A).*(p./t);
B

T = input('is maximum stress before shearing')
S=(T./B)