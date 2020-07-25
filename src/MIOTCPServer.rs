pub mod Server1
{
    use std::collections::HashMap;
    use std::io::{self, Read, Write};
    use std::str::from_utf8;

    use mio::{Events, Interest, Poll, Registry, Token};
    // You can run this example from the root of the mio repo:
// cargo run --example tcp_server --features="os-poll tcp"
    use mio::event::Event as Event1;
    use mio::net::{TcpListener, TcpStream};
    use quick_xml::Reader;
    use quick_xml::events::Event as Event2;


    // Setup some tokens to allow us to identify which event is for which socket.
    const SERVER: Token = Token(0);

    // Some data we'll send over the connection.
    const DATA: &[u8] = b"Hello world!\n";

    pub fn start() -> io::Result<()> {
        env_logger::init();

        // Create a poll instance.
        let mut poll = Poll::new()?;
        // Create storage for events.
        let mut events = Events::with_capacity(128);

        // Setup the TCP server socket.
        let addr = "192.168.1.4:60002".parse().unwrap();
        let mut server = TcpListener::bind(addr)?;

        // Register the server with poll we can receive events for it.
        poll.registry()
            .register(&mut server, SERVER, Interest::READABLE)?;

        // Map of `Token` -> `TcpStream`.
        let mut connections = HashMap::new();
        // Unique token for each incoming connection.
        let mut unique_token = Token(SERVER.0 + 1);

        println!("You can connect to the server using `nc`:");
        println!(" $ nc 192.168.1.4:60002");
        println!("You'll see our welcome message and anything you type we'll be printed here.");

        loop {
            poll.poll(&mut events, None)?;

            for event in events.iter() {
                match event.token() {
                    SERVER => loop {
                        // Received an event for the TCP server socket, which
                        // indicates we can accept an connection.
                        let (mut connection, address) = match server.accept() {
                            Ok((connection, address)) => (connection, address),
                            Err(e) if e.kind() == io::ErrorKind::WouldBlock => {
                                // If we get a `WouldBlock` error we know our
                                // listener has no more incoming connections queued,
                                // so we can return to polling and wait for some
                                // more.
                                break;
                            }
                            Err(e) => {
                                // If it was any other kind of error, something went
                                // wrong and we terminate with an error.
                                return Err(e);
                            }
                        };

                        println!("Accepted connection from: {}", address);

                        let token = next(&mut unique_token);
                        poll.registry().register(
                            &mut connection,
                            token,
                            Interest::READABLE.add(Interest::WRITABLE),
                        )?;

                        connections.insert(token, connection);
                    },
                    token => {
                        // Maybe received an event for a TCP connection.
                        let done = if let Some(connection) = connections.get_mut(&token) {
                            handle_connection_event(poll.registry(), connection, event)?
                        } else {
                            // Sporadic events happen, we can safely ignore them.
                            false
                        };
                        if done {
                            connections.remove(&token);
                        }
                    }
                }
            }
        }
    }

    fn next(current: &mut Token) -> Token {
        let next = current.0;
        current.0 += 1;
        Token(next)
    }

    /// Returns `true` if the connection is done.
    fn handle_connection_event(
        registry: &Registry,
        connection: &mut TcpStream,
        event: &Event1,
    ) -> io::Result<bool> {
        if event.is_writable() {
            // We can (maybe) write to the connection.
            match connection.write(DATA) {
                // We want to write the entire `DATA` buffer in a single go. If we
                // write less we'll return a short write error (same as
                // `io::Write::write_all` does).
                Ok(n) if n < DATA.len() => return Err(io::ErrorKind::WriteZero.into()),
                Ok(_) => {
                    // After we've written something we'll reregister the connection
                    // to only respond to readable events.
                    registry.reregister(connection, event.token(), Interest::READABLE)?
                }
                // Would block "errors" are the OS's way of saying that the
                // connection is not actually ready to perform this I/O operation.
                Err(ref err) if would_block(err) => {}
                // Got interrupted (how rude!), we'll try again.
                Err(ref err) if interrupted(err) => {
                    return handle_connection_event(registry, connection, event)
                }
                // Other errors we'll consider fatal.
                Err(err) => return Err(err),
            }
        }

        if event.is_readable() {
            let mut connection_closed = false;
            let mut received_data = Vec::with_capacity(4096);
            // We can (maybe) read from the connection.
            loop {
                let mut buf = [0; 256];
                match connection.read(&mut buf) {
                    Ok(0) => {
                        // Reading 0 bytes means the other side has closed the
                        // connection or is done writing, then so are we.
                        connection_closed = true;
                        break;
                    }
                    Ok(n) => received_data.extend_from_slice(&buf[..n]),
                    // Would block "errors" are the OS's way of saying that the
                    // connection is not actually ready to perform this I/O operation.
                    Err(ref err) if would_block(err) => break,
                    Err(ref err) if interrupted(err) => continue,
                    // Other errors we'll consider fatal.
                    Err(err) => return Err(err),
                }
            }

            if let Ok(str_buf) = from_utf8(&received_data) {
                let start=str_buf.find('<');
                if start.is_some()
                {
                    let start_index = start.unwrap();
                    let mut reader = Reader::from_str(&str_buf[start_index..]);
                    reader.trim_text(true);

                    let mut count = 0;
                    let mut txt = Vec::new();
                    let mut buf = Vec::new();

                    // The `Reader` does not implement `Iterator` because it outputs borrowed data (`Cow`s)
                    loop {
                        match reader.read_event(&mut buf) {
                            // for triggering namespaced events, use this instead:
                            // match reader.read_namespaced_event(&mut buf) {
                            Ok(Event2::Start(ref e)) => {
                                // for namespaced:
                                // Ok((ref namespace_value, Event::Start(ref e)))
                                match e.name() {
                                    b"tag1" => println!("attributes values: {:?}",
                                                        e.attributes().map(|a| a.unwrap().value)
                                                            .collect::<Vec<_>>()),
                                    b"tag2" => count += 1,
                                    _ => (),
                                }
                            },
                            // unescape and decode the text event using the reader encoding
                            Ok(Event2::Text(e)) => txt.push(e.unescape_and_decode(&reader).unwrap()),
                            Ok(Event2::Eof) => break, // exits the loop when reaching end of file
                            Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                            _ => (), // There are several other `Event`s we do not consider here
                        }

                        // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
                        buf.clear();
                    }

                    println!("Received data: {}", str_buf.trim_end());
                    for valss in txt {
                        println!("Received data2: {}", valss);
                    }
                }

            } else {
                println!("Received (none UTF-8) data: {:?}", &received_data);
            }

            if connection_closed {
                println!("Connection closed");
                return Ok(true);
            }
        }

        Ok(false)
    }

    fn would_block(err: &io::Error) -> bool {
        err.kind() == io::ErrorKind::WouldBlock
    }

    fn interrupted(err: &io::Error) -> bool {
        err.kind() == io::ErrorKind::Interrupted
    }
}