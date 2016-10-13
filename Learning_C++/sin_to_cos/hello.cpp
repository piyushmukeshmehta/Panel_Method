// Comment: prints "Hello world!" and an OS-independent newline
  #include <string>    // Defines type std::string
  #include <iostream>  // Defines global object std::cout
  using namespace std; // Allow std:: to be dropped
  int main() {         // Execution starts here
    string s="Hello world!\n"; // Declares object s of type string
    cout << s;         // An expression as a statement, << is the output operator
    return 0;          // Execution ends here
  }