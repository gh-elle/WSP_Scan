function bw = binarize_card(R, G, B, card_mask, threshold)
%BINARIZE_CARD  Binarise a WSP card by thresholding the red channel.
%
%   Stain pixels on WSP cards appear dark (low red reflectance) against the
%   bright yellow background.  The binary stain mask is:
%       bw = (R < threshold)  AND  card_mask
%
%   Parameters
%   ----------
%   R, G, B     : uint8 channel matrices for the full scan image
%   card_mask   : logical mask isolating the current card (size H × W)
%   threshold   : red-channel intensity cutoff (recommended: 120)
%
%   Returns
%   -------
%   bw : logical stain mask (size H × W; true = stain pixel)

bw = (R < threshold) & card_mask;

return
