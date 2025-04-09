% topoplotMEA - An EEGLAB function-like topography plotting function for
%               a NeuroNexus Mouse EEG 30-ch system.
% Usage:
%        >>  topoplotMEA(datavector, EEG.chanlocs, minMax); 
%
% Required Inputs:
%   datavector   - Single vector of channel values.
%   EEG.chanlocs - EEGLAB chanlocs structure. 
%   minMax       - Color scale range. If empty, use the data range.

% Copyright (C) Makoto Miyakoshi, Cincinnati Children's Hospital Medical Center
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% History
% 04/08/2025 Makoto. Created.

function targetAxes = topoplotMEA(data, chanlocs, minMax)

% Convert to double.
data = double(data);

% Define the xy coordinates.
x = [ 4.05  2.24  1.00  4.13  2.88  1.13  4.05  2.88  1.12  3.50  2.12 1.93 0.55 1.50 0.50 -0.50 -1.50 -0.55 -1.93 -2.12 -3.50 -1.12 -2.88 -4.05 -1.13 -2.88 -4.13 -1.00 -2.24 -4.05]; % 1x30 vector
y = [-4.14 -4.14 -4.14 -3.04 -3.04 -3.04 -1.96 -1.96 -1.96 -0.48 -0.48 1.04 1.04 2.30 2.30  2.30  2.30  1.04  1.04 -0.48 -0.48 -1.96 -1.96 -1.96 -3.04 -3.04 -3.04 -4.14 -4.14 -4.14]; % 1x30 vector

% Define interpolation grid (0.01 mm resolution)
[xq, yq] = meshgrid(-4.05:0.01:4.05, -4.14:0.01:2.30);

% Obtain sorting index.
chanIdx = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), {chanlocs.labels}');
[~, sortingIdx] = sort(chanIdx,'ascend');

% Interpolate signal values using cubic interpolation.
vq = griddata(x, y, data(sortingIdx), xq, yq, 'cubic');

% Determine the color scale.
if isempty(minMax)
    absMax = max(abs(vq(:)));
    minMax = [absMax*-1 absMax];
end

% Plot.
imageHandle = imagesc([-4.05 4.05], [-4.14 2.30], vq, minMax); 
set(gca, 'YDir', 'normal');  % ensure y increases upward
colormapMatrix = colormap('jet');
%colorbar;
%title('30-Ch MEA Topography');
xlabel('X (mm) to bregma');
ylabel('Y (mm) to bregma');
axis equal tight;

hold on

% Create a mask image (1 where NaN, 0 elsewhere).
nan_mask = isnan(vq);

% Create an RGB image to display the NaN overlay color.
overlay_color = [1 1 1];
nan_rgb = zeros([size(vq), 3]);  % preallocate RGB image.
for k = 1:3
    channel = overlay_color(k) * nan_mask;
    nan_rgb(:,:,k) = channel;
end

% Overlay it with full opacity where NaN, zero elsewhere.
h_overlay = imagesc([-4.05 4.05], [-4.14 2.30], nan_rgb);
set(h_overlay, 'AlphaData', nan_mask * 1);  % 1 = fully opaque.

% Add black disks of diameter 0.5 mm.
disk_diameter = 0.5; % mm

% Get CData (matrix of data values)
cdata = get(imageHandle, 'CData');

% Get colormap and clim (color limits)
clim = caxis; % [cmin, cmax]

% Normalize data to index into colormap
cmin = clim(1);
cmax = clim(2);
idx = round( (cdata - cmin) / (cmax - cmin) * (size(colormapMatrix,1)-1) ) + 1;

% Clamp indices to [1, size(cmap,1)]
idx = max(1, min(size(colormapMatrix,1), idx));

% Convert index matrix to RGB image
rgbImage = ind2rgb(idx, colormapMatrix);

for chIdx = 1:length(data)
    [~, rowIdx] = min(abs(yq(:,1) - y(sortingIdx(chIdx))));   % closest row index (y)
    [~, colIdx] = min(abs(xq(1,:) - x(sortingIdx(chIdx))));   % closest col index (x)
    
    % Extract RGB values at this coordinate
    if rowIdx > 1
        chRgb = squeeze(rgbImage(rowIdx-1, colIdx, :))'; % -1 was found heuristically. (04/08/2025 Makoto)
    else
        chRgb = squeeze(rgbImage(rowIdx, colIdx, :))';
    end

    rectangle('Position',[x(sortingIdx(chIdx))-disk_diameter/2, y(sortingIdx(chIdx))-disk_diameter/2, disk_diameter, disk_diameter],...
              'Curvature',[1,1],...
              'FaceColor', chRgb,...
              'EdgeColor', [0 0 0]);
end
hold off