startMessage <- function() {
  
  msg <- c(paste0(
  "     ______             __  __        _____ __  _____ __  ______
    / ____/         <_>/ / / /       / __  __ \\/ __  __ \\/ __  /
   / /_ __  ____    __/ /_/ /___  __/ / / / / / / / / / / /_/ /___________
  / ___/ _\\/ __ \\  / / /_   _/ / / / / / / / / / / / / / ____/ <_>_/ __  / 
 / /  / / / /_/ / / / / / / / /_/ / / / / / / / / / / / /   / /___/ / / /
/_/  /_/  \\___/\\\\/_/_/ /_/ /___  /_/  \\/  \\/_/  \\/  \\/_/   /_____/_/ /_/
                            _ / /   
                           \\___/ 
  Package version: ", packageVersion("frailtyMMpen")))
  
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  msg <- startMessage()
  if(!interactive())
    msg[1] <- paste0("Package frailtyMMpen version: ", packageVersion("frailtyMMpen"))
  packageStartupMessage(msg)      
  invisible()
}
