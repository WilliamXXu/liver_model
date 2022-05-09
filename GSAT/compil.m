% %rpi=raspi('192.168.1.153','pi','zx');
% led = rpi.AvailableLEDs{1};
% for i = 1:10
%     writeLED(rpi, led, 0);
%     pause(0.5);
%     writeLED(rpi, led, 1);
%     pause(0.5);
% end
% 
% system(rpi, 'ls -al /home/pi')
cfg = coder.config('lib','ecoder',false);
cfg.HardwareImplementation.ProdHWDeviceType = 'ARM Compatible->ARM Cortex';
cfg.Hardware = coder.hardware('Raspberry Pi'); 
cfg.GenCodeOnly = true;

codegen -config cfg sens -report

myBuildInfoFile = 'codegen/lib/sens/buildInfo.mat';
load(myBuildInfoFile);
packNGo(buildInfo);
% movefile ./codegen/lib/dashboard/examples/main.h
% movefile ./codegen/lib/dashboard/examples/main.c