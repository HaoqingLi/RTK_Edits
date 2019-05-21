classdef (ConstructOnLoad) AISEventDataPrediction < event.EventData
   properties
      dt
   end
   
   methods
      function data = AISEventDataPrediction(dt)
         data.dt = dt;
      end
   end
end