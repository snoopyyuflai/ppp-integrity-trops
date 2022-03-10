function obs = checkMeas(obj,obs0)
% This method is just a pass-through to the actual navigation filter

obs = obj.filter.checkMeas(obs0);

end