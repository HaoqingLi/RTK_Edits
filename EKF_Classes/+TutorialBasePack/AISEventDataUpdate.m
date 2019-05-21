classdef (ConstructOnLoad) AISEventDataUpdate < event.EventData
   properties
      obs
      dt
   end
   
   methods
      function data = AISEventDataUpdate(dt,obs)
         data.obs = obs;
         data.dt = dt;
      end
   end
end