// Dear emacs, this is -*- c++ -*-

#include "include/ObjectHandler.h"

#include <iostream>

ObjectHandler* ObjectHandler::m_instance = NULL;

ObjectHandler* ObjectHandler::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new ObjectHandler();
  }

  return m_instance;    

}

ObjectHandler::ObjectHandler() : m_logger( "ObjectHandler" )
{
  // constructor: initialise all variables
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;
  m_bcc = NULL;
  m_lumi = NULL;
}

ObjectHandler::~ObjectHandler()
{
  // default destructor

}

void ObjectHandler::SetBaseCycleContainer(BaseCycleContainer* bcc)
{
  // set the internal pointer to the container with all objects
  m_bcc = bcc;
  m_logger << DEBUG << "Pointer to BaseCycleContainer set." << SLogger::endmsg;
}

BaseCycleContainer* ObjectHandler::GetBaseCycleContainer()
{
  // return the pointer to the container with all objects
  if (!m_bcc){
    m_logger << WARNING << "Pointer to BaseCycleContainer is NULL." << SLogger::endmsg;
  }
  return m_bcc;
}

void ObjectHandler::SetLumiHandler(LuminosityHandler* lh)
{
  // set the internal pointer to the container with all objects
  m_lumi = lh;
  m_logger << DEBUG << "Pointer to LumiHandler set." << SLogger::endmsg;
}

LuminosityHandler* ObjectHandler::GetLumiHandler()
{
  // return the pointer to the container with all objects
  if (!m_lumi){
    m_logger << WARNING << "Pointer to LumiHandler is NULL." << SLogger::endmsg;
  }
  return m_lumi;
}

