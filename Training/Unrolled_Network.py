import torch
from DC import DC_2d
from ResNet import ResNet

class UnrolledNet(torch.nn.Module):
    
    def __init__(self, mu = 0.0, Unrolls = 10, device="cuda:3", echoes = 5):
        super().__init__()
        self.Network         = ResNet(in_Chn = echoes*2, inner_Chn = 64, out_Chn = echoes*2)  
        self.DataConsistency = DC_2d(mu, iters=15, echos = echoes)
        self.Unrolls         = Unrolls
        self.echos           = echoes
        self.dev             = device


    def complex2real(self,x): # (1,6,120,120) --> (1,6,2,120,120)
        x = x.unsqueeze(2)
        return torch.concat((torch.real(x),torch.imag(x)), dim=2)

    def real2complex(self,x): # (1,6,2,120,120) --> (1,6,120,120)
        return x[: , : , 0 , : , :] + 1j * x[: , : , 1 , : , :]

    # # Without CheckPointing
    def forward(self, zf , coil_map, M):
        cg    = DC_2d(mu=0, iters=3, echos=self.echos, flag_learnable=False)
        recon = cg(zf, zf, coil_map, M)
        
        for OuterIter in range(self.Unrolls):
            recon = self.complex2real(recon)       # [B, E, 2, H, W]
            recon = self.Network(recon.float())    # [B, 2E, H, W]
            recon = self.real2complex(recon)       # [B, E, H, W]
            recon = self.DataConsistency(recon, zf, coil_map, M)
        return recon
    
class UnrolledNet_PE(torch.nn.Module):
    
    def __init__(self, mu = 0.0, Unrolls = 10, device="cuda:3", echoes = 5):
        super().__init__()
        self.Network         = ResNet(in_Chn = echoes*2, inner_Chn = 64, out_Chn = echoes*2)  
        self.DataConsistency = DC_2d(mu, iters=15, echos = echoes)
        self.Unrolls         = Unrolls
        self.echos           = echoes
        self.dev             = device


    def complex2real(self,x): # (1,6,120,120) --> (1,6,2,120,120)
        x = x.unsqueeze(2)
        return torch.concat((torch.real(x),torch.imag(x)), dim=2)

    def real2complex(self,x): # (1,6,2,120,120) --> (1,6,120,120)
        return x[: , : , 0 , : , :] + 1j * x[: , : , 1 , : , :]

    # # Without CheckPointing
    def forward(self, zf , coil_map, M, p_prime):
        cg    = DC_2d(mu=0, iters=3, echos=self.echos, flag_learnable=False)
        recon = cg(zf, zf, coil_map, M)
        
        for OuterIter in range(self.Unrolls):
            recon = recon * torch.conj(p_prime)    # Phase Corection
            recon = self.complex2real(recon)       # [B, E, 2, H, W]
            recon = self.Network(recon.float())    # [B, 2E, H, W]
            recon = self.real2complex(recon)       # [B, E, H, W]
            recon = recon * p_prime                # Phase Corection
            recon = self.DataConsistency(recon, zf, coil_map, M)
        return recon
