function plotTvRecovery(problem, results, savePrefix)
%PLOTTVRECOVERY  精确间距（15pt）的图像对比图，支持 LaTeX 下标，PSNR 括号显示。
%
%   布局参数（单位：磅 points）：
%       每子图轴尺寸 180×160 pt，水平间距 15 pt，垂直间距 15 pt。
%       总子图 ≤3 时单行，否则双行。
%   字体：Times New Roman 12pt；标题使用 LaTeX 解释器。
%   算法名中仅 _{spp} 被替换为 $_{\rm spp}$ 以显示下标，其余名称不变。
%   噪声标题为 “Noisy (PSNR dB)”，算法标题为 “算法名 (PSNR dB)”。
%
%   输入：
%       problem    - 含 xtrue, m, n, c 的结构体
%       results    - 元胞数组，每个元素含 .x（图像）和 .name（算法名）
%       savePrefix - 保存路径前缀（可选）
%   输出：300 DPI 的 PNG 文件，并关闭图形窗口。

% ==================== 默认路径 ====================
if nargin < 3 || isempty(savePrefix)
    savePrefix = defaultSavePrefix('tv_recovery');
end

% ==================== 图像数据 ====================
origImg  = reshape(problem.xtrue, problem.m, problem.n);
noisyImg = reshape(problem.c,     problem.m, problem.n);
numAlgs   = numel(results);
numImages = 2 + numAlgs;          % 原始 + 噪声 + 各算法

% ==================== 自适应行数 ====================
if numImages <= 3
    numRows = 1;
else
    numRows = 2;
end
numCols = ceil(numImages / numRows);

% ==================== 布局参数（全部以 points 为单位） ====================
axesWidth  = 180;   % 每个子图轴的宽度
axesHeight = 160;   % 每个子图轴的高度
hSpacing   = 15;    % 水平间距
vSpacing   = 18;    % 垂直间距（上下排之间的距离）
leftMargin = 30;    % 左边距
bottomMargin = 40;  % 底边距
topMargin = 20;     % 顶边距（为标题预留空间）
rightMargin = 20;   % 右边距

% 计算整个图形窗口的总尺寸（points）
totalWidth  = leftMargin + numCols*axesWidth  + (numCols-1)*hSpacing + rightMargin;
totalHeight = bottomMargin + numRows*axesHeight + (numRows-1)*vSpacing + topMargin;

% ==================== 创建图形窗口（单位：points） ====================
fig = figure('Color', 'w', 'Name', 'TV recovery comparison', ...
             'Visible', 'off', 'Units', 'points', ...
             'Position', [100, 100, totalWidth, totalHeight]);

% 全局字体设置
set(fig, 'DefaultTextFontName', 'Times New Roman');
set(fig, 'DefaultAxesFontName', 'Times New Roman');
set(fig, 'DefaultTextFontSize', 12);
set(fig, 'DefaultAxesFontSize', 12);

% ==================== 绘制子图 ====================
imgIdx = 0;   % 从 0 开始计数，方便行列计算

% ---- 1. 原始图像 ----
imgIdx = imgIdx + 1;
plotOneImage(fig, origImg, 'Original', imgIdx, numCols, numRows, ...
             axesWidth, axesHeight, hSpacing, vSpacing, ...
             leftMargin, bottomMargin, topMargin, totalHeight);

% ---- 2. 噪声图像 ----
imgIdx = imgIdx + 1;
psnrNoisy = Functions.psnrValue(problem, problem.c);
noiseTitle = sprintf('Noisy (%.2f dB)', psnrNoisy);
plotOneImage(fig, noisyImg, noiseTitle, imgIdx, numCols, numRows, ...
             axesWidth, axesHeight, hSpacing, vSpacing, ...
             leftMargin, bottomMargin, topMargin, totalHeight);

% ---- 3 ~ N. 各算法结果 ----
for i = 1:numAlgs
    imgIdx = imgIdx + 1;
    recImg = reshape(results{i}.x, problem.m, problem.n);
    psnrVal = Functions.psnrValue(problem, results{i}.x);

    % 处理算法名：仅替换 _{spp} 为 LaTeX 下标
    algName = results{i}.name;
    if contains(algName, '_{spp}')
        algNameDisplay = strrep(algName, '_{spp}', '$_{\rm spp}$');
    else
        algNameDisplay = algName;
    end
    algTitle = [algNameDisplay, sprintf(' (%.2f dB)', psnrVal)];

    plotOneImage(fig, recImg, algTitle, imgIdx, numCols, numRows, ...
                 axesWidth, axesHeight, hSpacing, vSpacing, ...
                 leftMargin, bottomMargin, topMargin, totalHeight);
end

% ==================== 保存与清理 ====================
outDir = fileparts(savePrefix);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
pngFile = [savePrefix '_recovery_comparison.png'];
set(fig, 'Renderer', 'painters');
exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
end

% ==================== 局部函数：绘制单个子图 ====================
function plotOneImage(fig, img, titleStr, imgIdx, numCols, numRows, ...
                      axesWidth, axesHeight, hSpacing, vSpacing, ...
                      leftMargin, bottomMargin, topMargin, totalHeight)
% plotOneImage 在指定网格位置显示图像并添加标题。
%
%   网格顺序：从左到右、从上到下。
%   imgIdx 从 1 开始。
%   行列计算：col = mod(imgIdx-1, numCols) + 1;
%            row = ceil(imgIdx / numCols);   (1 为最上方行)
%
%   坐标计算（单位：points）：
%       left   = leftMargin + (col-1)*(axesWidth + hSpacing)
%       bottom = totalHeight - topMargin - row*axesHeight - (row-1)*vSpacing

    col = mod(imgIdx-1, numCols) + 1;          % 列号（从左开始）
    row = ceil(imgIdx / numCols);              % 行号（从上到下：1 为第一行）

    % 轴左下角位置
    axLeft   = leftMargin + (col-1) * (axesWidth + hSpacing);
    axBottom = totalHeight - topMargin - row*axesHeight - (row-1)*vSpacing;

    % 创建坐标轴
    ax = axes('Units', 'points', ...
              'Position', [axLeft, axBottom, axesWidth, axesHeight], ...
              'Parent', fig);
    imshow(img, [], 'Parent', ax);
    axis(ax, 'image', 'off');

    % 标题（LaTeX 解释器，Times 12pt）
    title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 12);
end

% ==================== 局部函数：默认保存前缀 ====================
function savePrefix = defaultSavePrefix(baseName)
    toolDir = fileparts(mfilename('fullpath'));
    projectDir = fileparts(toolDir);
    savePrefix = fullfile(projectDir, 'Results', baseName);
end