mod MIOTCPServer;
mod TKIOTCPServer;
mod UDPServer;
pub use crate::MIOTCPServer::Server1;
pub use crate::UDPServer::Server2;
pub use crate::TKIOTCPServer::Server3;

fn main() {
    let i =1;
    if i==1
    {
        Server1::start();
    }
    else if i==2
    {
        Server2::start();
    }
    else if i==3
    {
        Server3::start();
    }

}