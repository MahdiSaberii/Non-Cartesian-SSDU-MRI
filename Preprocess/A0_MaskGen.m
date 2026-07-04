clear all;close all;clc;
TotalNumberofPoints = 3592;
Mask_n              = 7;
ssdu_ratio          = 0.4; 
save_path           = './Masks/';

if ~exist(save_path, 'dir')
    mkdir(save_path)
end

for ii = 1: Mask_n
    count = round(TotalNumberofPoints*ssdu_ratio);
    validation = zeros(1,TotalNumberofPoints);
while count > 0
    a = round(rand(1, 1)*TotalNumberofPoints);
    if a == 0
        a = TotalNumberofPoints;
    end
    if validation(a) == 0 && a > 32
        validation(a) = 1;
        count = count - 1;
    end
end

training = ~validation;
training = double(training);

file_name = sprintf('%smask_%d.mat', save_path, ii);
save(file_name, "training", "validation")

end
