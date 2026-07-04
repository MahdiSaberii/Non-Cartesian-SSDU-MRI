import sys
import torch

class DC_2d(torch.nn.Module):
 
    def __init__(self, mu=0.0, iters = 15, echos = 5, flag_learnable = True):
        super().__init__()
        self.mu    = torch.nn.Parameter(torch.ones([echos],dtype = torch.float32)*mu, requires_grad=flag_learnable)
        self.iters = iters
        self.echos = echos
        
    def M_operator (self , image, coil , M):
        coil_images = image * coil
        size_x = coil_images.shape[1] // 2
        size_y = coil_images.shape[2] // 2
        
        a = torch.nn.functional.pad(coil_images,[size_x,size_x,size_y,size_y])
        e = torch.fft.fftn(torch.fft.fftshift(a, dim=[-2,-1]),dim=[-2,-1], norm='ortho')
        f = e * M
        g = torch.fft.ifftshift(torch.fft.ifftn(f ,dim=[-2,-1] ,norm='ortho'), dim=[-2,-1])
        h = g[:,size_x:size_x*3,size_y:size_y*3]
        out = torch.sum(h*torch.conj(coil),axis=0)
        
        return out
    
    def forward(self, cnn_out , zerofilled ,coil, M): # M: [1,n_echo,120,120]

        epsilon          = torch.tensor(sys.float_info.epsilon)
        all_out          = torch.zeros_like(zerofilled)  # [1,n_echo,120,120]
        all_b_iterations = torch.zeros((self.iters, *zerofilled.shape), dtype=torch.complex64)  # [iters, 1, n_echo, 120, 120]

        for echo_idx in range(self.echos):
            p = zerofilled[:, echo_idx, :, :] + self.mu[echo_idx] * cnn_out[:, echo_idx, :, :]
            b = torch.zeros_like(p)
            r = torch.clone(p)
            
            for i in range(self.iters):
                
                q     = self.M_operator(p, coil, M[:,echo_idx, :, :]) + self.mu[echo_idx] * p
                rsold = (r * torch.conj(r))
                alpha = (torch.sum(rsold) / torch.sum(q * torch.conj(p) + epsilon))
                b     = b + alpha * p
                r     = r - alpha * q
                rsnew = (r * torch.conj(r))
                p     = r + (torch.sum(rsnew) / torch.sum(rsold + epsilon)) * p
                
                all_b_iterations[i, :, echo_idx, :, :] = b
            all_out[:, echo_idx, : , :] = b
        return all_out
