load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/15cg_7masks/testresultsb6_x10_ssdu1.mat'

recon1 = squeeze(recon);

load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/15cg_7masks/testresultsb6_x10_ssdu2.mat'

recon2 = squeeze(recon);

load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/15cg_7masks/testresultsb6_x10_ssdu3.mat'

recon3 = squeeze(recon);

%%
figure()
for i = 1:6
    subplot(4,6,i+6)
    imshow(squeeze(abs(recon1(i,:,:))), [])
    
    subplot(4,6,i+12)
    imshow(squeeze(abs(recon2(i,:,:))), [])
    
    subplot(4,6,i+18)
    imshow(squeeze(abs(recon3(i,:,:))), [])
end

%%
load '/home/naxos2-raid7/hongygu/non_cartesian/spiral_fmri/DCF_modified_results/my_training/ssdu_0_5/new_run/fixed_data/fixed/15cg_7masks/reg10.mat'
e = 1;
i = 1;
figure()
subplot(4,10,i)
imshow(squeeze(abs(ref_2(i,:,:))), [])
    
subplot(4,10,i+10)
imshow(squeeze(abs(atb_2(i,:,:))), [])
    
subplot(4,10,i+20)
imshow(squeeze(abs(recon1(i,:,:))), [])
    
subplot(4,10,i+30)
imshow(squeeze(abs(recon1(i,:,:))), [])