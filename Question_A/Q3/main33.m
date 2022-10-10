clear all;clc;close all;

Ts=1.25e-7;
T=3.2E-5;

N=T/Ts;
Nf=32;
L=0.0815;
gamma=78.986e12;
Na=86;
f0=78.8e9;
c=299792458;
lamda=c/f0;

MaxARange=55;
MinARange=8;
MaxRRange=15;
MinRRange=1;
Resl=0.05;
AntDis=L/(Na-1);

load('../data/data_q3.mat');