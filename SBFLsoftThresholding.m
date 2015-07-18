function out = SBFLsoftThresholding(in, t)
    out = in;
    idx1 = out < -t;
    idx2 = out > t;
    out(~idx1 & ~idx2) = 0;
    out(idx1) = out(idx1) + t;
    out(idx2) = out(idx2) - t;
end