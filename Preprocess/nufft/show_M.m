load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/all_ones_mask.mat'

load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/0_4/M_0_4.mat'

M_train = squeeze(M(1,:,:,:));

M_test = squeeze(M(6,:,:,:));

%%
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu/training_slices_simple/Fessler_ssdu_M_simple_testing_x10'

M_10 = M;

load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/M_ref.mat'

M_full = M;

for i = 1: 6
    M_train(i,:,:) = squeeze(M_train(i,:,:)) + mask_M;
    M_test(i,:,:) = squeeze(M_test(i,:,:)) + mask_M;
    M_10(i,:,:) = squeeze(M_10(i,:,:)) + mask_M;
    M_full(i,:,:) = squeeze(M_full(i,:,:)) + mask_M;
end

p1 = [squeeze(M_full(1,:,:)), squeeze(M_full(2,:,:)), squeeze(M_full(3,:,:)),squeeze(M_full(4,:,:)),squeeze(M_full(5,:,:)),squeeze(M_full(6,:,:))];
p2 = [squeeze(M_10(1,:,:)), squeeze(M_10(2,:,:)), squeeze(M_10(3,:,:)),squeeze(M_10(4,:,:)),squeeze(M_10(5,:,:)),squeeze(M_10(6,:,:))];
p3 = [squeeze(M_train(1,:,:)), squeeze(M_train(2,:,:)), squeeze(M_train(3,:,:)),squeeze(M_train(4,:,:)),squeeze(M_train(5,:,:)),squeeze(M_train(6,:,:))];
p4 = [squeeze(M_test(1,:,:)), squeeze(M_test(2,:,:)), squeeze(M_test(3,:,:)),squeeze(M_test(4,:,:)),squeeze(M_test(5,:,:)),squeeze(M_test(6,:,:))];

%%
p = [p1;p2;p3;p4];

figure()
imshow(abs(p), [0 1])

%%
imwrite(abs(squeeze(M_10(1,:,:))), '/home/naxos2-raid7/hongygu/non_cartesian/M_example.png')

%%




