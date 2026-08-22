function bank = make_sharing_template_bank_v06(blocks,innerUniverse,termsPerTemplate)
%MAKE_SHARING_TEMPLATE_BANK_V06 Deterministic equal-cardinality block masks.
%
% Row 1 is deliberately highly structured: the first K low exponents.
% Rows 2..blocks are deterministic pseudo-random-looking K-of-U masks made
% from a fixed sine hash.  No RNG state or toolbox is required, making the
% experiment reproducible across repeated runs.
if blocks < 1 || blocks ~= floor(blocks)
    error('BMD:V06Bank','blocks must be a positive integer.');
end
if innerUniverse < 2 || innerUniverse ~= floor(innerUniverse)
    error('BMD:V06Bank','innerUniverse must be an integer >=2.');
end
if termsPerTemplate < 1 || termsPerTemplate >= innerUniverse || termsPerTemplate ~= floor(termsPerTemplate)
    error('BMD:V06Bank','termsPerTemplate must be an integer in 1..innerUniverse-1.');
end

bank = false(blocks,innerUniverse);
bank(1,1:termsPerTemplate) = true;
i = 1:innerUniverse;
for t = 2:blocks
    scores = mod(sin(i*12.9898 + t*78.233)*43758.5453123,1);
    [~,ord] = sort(scores,'ascend');
    bank(t,ord(1:termsPerTemplate)) = true;
end

if any(sum(bank,2) ~= termsPerTemplate)
    error('BMD:V06Bank','Template cardinality invariant failed.');
end
if size(unique(bank,'rows'),1) ~= blocks
    error('BMD:V06Bank','Template generator produced duplicate rows.');
end
end
