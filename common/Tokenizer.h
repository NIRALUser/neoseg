
////////////////////////////////////////////////////////////////////////////////
//
// Simple string tokenizer
//
////////////////////////////////////////////////////////////////////////////////


#ifndef _Tokenizer_h
#define _Tokenizer_h

class Tokenizer {

private:
  char* m_str;
  int m_str_len;
  char* m_curr_token;
  char* m_delims;
  int m_delims_len;
  int m_offt;

public:
  Tokenizer();
  ~Tokenizer();

  void SetDelimiters(const char* s);
  void SetString(const char* s);

  // Returns NULL when no more tokens are available
  char* GetNextToken();

  // Returns the rest of the untokenized string
  char* GetRemainder();

};

#endif
