from matplotlib.axes._axes import Axes
import numpy as np

def get_extent(img:np.ndarray,ax:Axes,L:float,pc:tuple[float,float],W:float,H:float):
    """
    获取显示图像的参数
    img
    """
    h,w,_=img.shape
    # W,H=8,6
    l=ax.figure.subplotpars.left
    r=ax.figure.subplotpars.right
    t=ax.figure.subplotpars.top
    b=ax.figure.subplotpars.bottom
    x0,x1=ax.get_xlim()
    y0,y1=ax.get_ylim()
    K=(w*H*(t-b)*(x1-x0))/(h*W*(r-l)*(y1-y0))
    xc,yc=pc #图片中心点坐标
    x0_,x1_,y0_,y1_=xc-L/2, xc+L/2, yc-L/(2*K), yc+L/(2*K)
    return x0_,x1_,y0_,y1_