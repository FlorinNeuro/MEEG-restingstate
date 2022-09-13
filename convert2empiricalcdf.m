function [cdf_values] = convert2empiricalcdf(data_vector,include0)

% This function computes the emperical cummulative density function for the
% data_vector. Need to specifey if 0 should be included or not.

% Inputs:
% pdf = probability density function of the data_vector (normalized/non-normalized)
% pdf = (number_of_samples X 1)
% data_vector = data associated with the pdf
% data_vector = (number_of_samples X 1)
% include 0 wether or not 0 should be considered for the empirical
% distribution (SVD approaches currently not included, ICA yes)

% Esther Florin 2022



if include0==0
    h=histogram(data_vector(data_vector>0),size(data_vector(data_vector>0), 1),'Normalization','probability' );
else
    h=histogram(data_vector,size(data_vector, 1),'Normalization','probability' );
end
% normalization
pdf = h.Values;
cdf_values = zeros(length(pdf),1);
if include0==0
    for dt = 1:length(data_vector)

        less_data_points = find(data_vector (data_vector>0)< data_vector(dt));

        cdf_values(dt,1) = sum(pdf(less_data_points));

    end
else
    for dt = 1:length(data_vector)
        less_data_points = find(data_vector < data_vector(dt));
        cdf_values(dt,1) = sum(pdf(less_data_points));
    end

end

% less_data_points = find(data_vector (data_vector>0)< data_point);

% Plotting calculated cdf
%     figure
%     set(gcf,'Position',[440   378   250   400])
%     [cdf_sorted,I_sorted] = sort(cdf_values);
%     data_vector_sorted = data_vector(I_sorted);
%     plot(data_vector_sorted,cdf_sorted,'b','LineWidth',4)
end