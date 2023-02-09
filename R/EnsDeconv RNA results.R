# EnsDeconv RNA results
load("D:/Manqi/abstract/0415/0415archive.RData")
test <- res_all_sca %>% filter(str_detect(Method,"EnsDeconv")& bulk == "ROS") %>% ungroup() %>% select(Subject,CellType,p_hat)
ros_5frac <- spread(test, key = CellType, value = p_hat)%>% remove_rownames %>% column_to_rownames(var="Subject")
ros_5frac <- as.matrix(ros_5frac)
saveRDS(ros_5frac,"EnsDeconv_ROS_result.rds")
