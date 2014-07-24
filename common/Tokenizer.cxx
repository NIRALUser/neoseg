
#include "Tokenizer.h"

#include <iostream>

#include <string.h>

Tokenizer::
Tokenizer()
{
  m_str = NULL;
  m_str_len = 0;
  m_curr_token = NULL;
  // Whitespace as default delims
  char* s = " \t\n";
  m_delims = new char[strlen(s)+1];
  strcpy(m_delims, s);
  m_delims_len = strlen(m_delims);
  m_offt = 0;
}

Tokenizer::
~Tokenizer()
{
  if (m_str != NULL)
    delete [] m_str;
  if (m_curr_token != NULL)
    delete [] m_curr_token;
  if (m_delims != NULL)
    delete [] m_delims;
}

void
Tokenizer::
SetDelimiters(const char* s)
{

  if (s == NULL)
  {
    std::cerr << "[Tokenizer::SetDelimiters] NULL argument" << std::endl;
    return;
  }

  if (m_delims != NULL)
    delete [] m_delims;

  m_delims_len = strlen(s);
  m_delims = new char[m_delims_len+1];
  strcpy(m_delims, s);

  // Reset current status
  m_offt = 0;
  if (m_curr_token != NULL)
    m_curr_token[0] = 0;

}

void
Tokenizer::
SetString(const char* s)
{

  if (s == NULL)
  {
    std::cerr << "[Tokenizer::SetString] NULL argument" << std::endl;
    return;
  }

  unsigned int n = strlen(s);

  if (m_str != NULL)
    delete [] m_str;
  if (m_curr_token != NULL)
    delete [] m_curr_token;

  m_str = new char[n+1]; 
  m_str_len = n;
  strcpy(m_str, s);

  m_curr_token = new char[n+1];
  m_curr_token[0] = 0;

  m_offt = 0;

}

char*
Tokenizer::
GetNextToken()
{

  if (m_str == NULL)
    return NULL;

  if (m_offt >= m_str_len)
    return NULL;

  m_curr_token[0] = 0;

  int traversed = 0;
  int tok_offt = 0;

  bool valid = false;

  for (int i = m_offt; i < m_str_len; i++)
  {

    bool isdelim = false;
    for (int j = 0; j < m_delims_len; j++)
    {
      if (m_str[i] == m_delims[j]) {
        isdelim = true;
        break;
      }
    }

    traversed++;

    // Already have a token and now seen a delimiter
    if (valid && isdelim)
      break;

    // Place the char in the token and mark token as valid (not empty)
    if (!isdelim) {
      valid = true;
      m_curr_token[tok_offt] = m_str[i];
      tok_offt++;
    }

  }

  m_offt += traversed;

  if (!valid)
    return NULL;

  // Set last byte to zero
  m_curr_token[tok_offt] = 0;

  return m_curr_token;

}

char*
Tokenizer::
GetRemainder()
{

  // TODO
  // Returns the remaining substring

  return NULL;

}

// Test
#if 0

void print_tokens(Tokenizer* tok)
{

  std::cout << "****************" << std::endl;

  char* t = tok->GetNextToken();
  while (t != NULL) {
    std::cout << t << "!" << std::endl;
    t = tok->GetNextToken();
  }
  std::cout << std::endl;

}

int
main(int argc, char** argv)
{

  char* s1 = "abc def ghi ";
  char* s2 = "abc =  314.2";
  char* s3 = "1.3 1.0 0.7";

  Tokenizer* tok = new Tokenizer();

  tok->SetString(s1);
  print_tokens(tok);

  tok->SetString(s2);
  print_tokens(tok);

  tok->SetDelimiters("\t =");

  tok->SetString(s1);
  print_tokens(tok);

  tok->SetString(s2);
  print_tokens(tok);

  tok->SetDelimiters(" \t");

  tok->SetString(s3);
  print_tokens(tok);

  return 0;

}

#endif
