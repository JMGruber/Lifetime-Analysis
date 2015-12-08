function [HistFiltered] = FilterHistograms(Histograms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    index=max(Histograms(2:end,:).');
    index=[0 index];
    %selection=find(index>15);
    HistFiltered=Histograms(index>15,:);




